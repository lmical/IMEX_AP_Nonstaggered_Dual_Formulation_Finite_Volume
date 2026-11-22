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



#ifdef CENTEREDPRIMITIVE
#ifdef IMEX
USE MOD_Equation           ,ONLY: Compute_min_p_max_ro
#endif
#endif

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D     ,ONLY: Overwrite_WC_from_U
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
#ifndef DONOTOVERWRITEPRIMITIVE
  CALL Overwrite_WC_from_U()
#endif
#endif

#ifdef CENTEREDPRIMITIVE
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
! #ifdef CENTEREDPRIMITIVE
!     PRINT*, "WC ", iVar, MAXVAL(WC(iVar,1:nElemsX,1:nElemsY)), MINVAL(WC(iVar,1:nElemsX,1:nElemsY))
! #endif
!   END DO
!   PRINT*
  !*=========================================================


  SELECT CASE(timescheme)
#if defined(IMEX)
#ifdef CENTEREDPRIMITIVE
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
#ifdef CENTEREDPRIMITIVE
    CASE(-42)
      CALL TimeDiscretizationByARS222_Crank_Nicolson_Conserved(t)
    CASE(-82)
      CALL IMEXDeC2_rho_p_stage_dependent_overwriting_stages(t)
#endif
#endif

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

#endif
    CASE default
      WRITE(*,*) "Time discretization not implemented"
      STOP
  END SELECT



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
#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: K0_WC
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
#ifdef CENTEREDPRIMITIVE
K0_WC(1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
#endif

!--------------------!
! Forward Euler      !
!--------------------!
tStage = t 
CALL FVTimeDerivative(tStage)

U(1:nVar,1:nElemsX,1:nElemsY) = K0(1:nVar,1:nElemsX,1:nElemsY) + Ut(1:nVar,1:nElemsX,1:nElemsY)*dt
#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = K0_WC(1:nVar,1:nElemsX,1:nElemsY) + WCt(1:nVar,1:nElemsX,1:nElemsY)*dt
#endif


!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByForwardEuler
!===============================================================================!
!
!
!
!===============================================================================!
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
#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: K0_WC
USE MOD_FiniteVolume2D_vars,ONLY: K1_WC
USE MOD_FiniteVolume2D_vars,ONLY: K2_WC
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
#ifdef CENTEREDPRIMITIVE
K0_WC(1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
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
#ifdef CENTEREDPRIMITIVE
K1_WC(1:nVar,1:nElemsX,1:nElemsY) = &
    K0_WC(1:nVar,1:nElemsX,1:nElemsY) &
  + WCt(1:nVar,1:nElemsX,1:nElemsY)*dt
WC(1:nVar,1:nElemsX,1:nElemsY)  = K1_WC(1:nVar,1:nElemsX,1:nElemsY)
#endif



!--------------------!
! Second Stage       !
!--------------------!
tStage = t + 1.0*dt
CALL FVTimeDerivative(tStage)
K2(1:nVar,1:nElemsX,1:nElemsY) = &
    0.5*K0(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.5*K1(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.5*Ut(1:nVar,1:nElemsX,1:nElemsY)*dt
U(1:nVar,1:nElemsX,1:nElemsY)  = K2(1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
K2_WC(1:nVar,1:nElemsX,1:nElemsY) = &
    0.5*K0_WC(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.5*K1_WC(1:nVar,1:nElemsX,1:nElemsY) &
  + 0.5*WCt(1:nVar,1:nElemsX,1:nElemsY)*dt
WC(1:nVar,1:nElemsX,1:nElemsY)  = K2_WC(1:nVar,1:nElemsX,1:nElemsY)
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
#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: K0_WC
USE MOD_FiniteVolume2D_vars,ONLY: K1_WC
USE MOD_FiniteVolume2D_vars,ONLY: K2_WC
USE MOD_FiniteVolume2D_vars,ONLY: K3_WC
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
#ifdef CENTEREDPRIMITIVE
K0_WC(1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
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
#ifdef CENTEREDPRIMITIVE
K1_WC(1:nVar,1:nElemsX,1:nElemsY)= WCt(1:nVar,1:nElemsX,1:nElemsY)
WC(1:nVar,1:nElemsX,1:nElemsY) = &
    K0_WC(1:nVar,1:nElemsX,1:nElemsY) &
  + K1_WC(1:nVar,1:nElemsX,1:nElemsY)*dt
#endif


!--------------------!
! Second Stage       !
!--------------------!
tStage = t + 1.0*dt
CALL FVTimeDerivative(tStage)
K2(1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)
U(1:nVar,1:nElemsX,1:nElemsY) = K0(1:nVar,1:nElemsX,1:nElemsY) + &
          0.25*K1(1:nVar,1:nElemsX,1:nElemsY)*dt+0.25*K2(1:nVar,1:nElemsX,1:nElemsY)*dt 
#ifdef CENTEREDPRIMITIVE
K2_WC(1:nVar,1:nElemsX,1:nElemsY) = WCt(1:nVar,1:nElemsX,1:nElemsY)
WC(1:nVar,1:nElemsX,1:nElemsY) = K0_WC(1:nVar,1:nElemsX,1:nElemsY) + &
          0.25*K1_WC(1:nVar,1:nElemsX,1:nElemsY)*dt+0.25*K2_WC(1:nVar,1:nElemsX,1:nElemsY)*dt 
#endif

!--------------------!
! Third Stage       !
!--------------------!
tStage = t + 0.5*dt
CALL FVTimeDerivative(tStage)
K3(1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)
U(1:nVar,1:nElemsX,1:nElemsY) = K0(1:nVar,1:nElemsX,1:nElemsY) + &
          1./6.*K1(1:nVar,1:nElemsX,1:nElemsY)*dt+&
          1./6.*K2(1:nVar,1:nElemsX,1:nElemsY)*dt+&
          2./3.*K3(1:nVar,1:nElemsX,1:nElemsY)*dt
#ifdef CENTEREDPRIMITIVE
K3_WC(1:nVar,1:nElemsX,1:nElemsY) = WCt(1:nVar,1:nElemsX,1:nElemsY)
WC(1:nVar,1:nElemsX,1:nElemsY) = K0_WC(1:nVar,1:nElemsX,1:nElemsY) + &
          1./6.*K1_WC(1:nVar,1:nElemsX,1:nElemsY)*dt+&
          1./6.*K2_WC(1:nVar,1:nElemsX,1:nElemsY)*dt+&
          2./3.*K3_WC(1:nVar,1:nElemsX,1:nElemsY)*dt
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


#ifdef CENTEREDPRIMITIVE
  Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
  Wp_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
#endif

END DO


CALL FVTimeDerivative(tStage(1))

DO ii = 1,MSteps
  FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)


#ifdef CENTEREDPRIMITIVE
  FWp_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WCt(1:nVar,1:nElemsX,1:nElemsY)
#endif

END DO

DO kk = 1,KCorr
  IF (kk == KCorr) THEN
    ii = MSteps
    Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(1,1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(1,1:nVar,1:nElemsX,1:nElemsY)
#endif

    DO jj = 1,MSteps
      Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
      Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,1:nElemsX,1:nElemsY)
#endif
    END DO
  ELSE
    DO ii = 2,MSteps
      Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(1,1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
      Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(1,1:nVar,1:nElemsX,1:nElemsY)
#endif

      DO jj = 1,MSteps
        Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
        Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,1:nElemsX,1:nElemsY)
#endif
      END DO
    END DO
  END IF

  IF (kk .NE. KCorr) THEN
    DO ii = 2,MSteps
      U(1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
      WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY)
#endif

      CALL FVTimeDerivative(tStage(ii))
      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WCt(1:nVar,1:nElemsX,1:nElemsY)
#endif

    END DO
  END IF

  Up=Ua

#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

END DO

U(1:nVar,1:nElemsX,1:nElemsY) = Ua(MSteps,1:nVar,1:nElemsX,1:nElemsY)
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


#ifdef CENTEREDPRIMITIVE
  Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
  Wp_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
#endif

END DO


CALL FVTimeDerivative(tStage(1))

DO ii = 1,MSteps
  FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)


#ifdef CENTEREDPRIMITIVE
  FWp_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WCt(1:nVar,1:nElemsX,1:nElemsY)
#endif

END DO

DO kk = 1,KCorr
  IF (kk == KCorr) THEN
    ii = MSteps
    Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(1,1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(1,1:nVar,1:nElemsX,1:nElemsY)
#endif

    DO jj = 1,MSteps
      Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
      Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,1:nElemsX,1:nElemsY)
#endif
    END DO
  ELSE
    DO ii = 2,MSteps
      Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(1,1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
      Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(1,1:nVar,1:nElemsX,1:nElemsY)
#endif

      DO jj = 1,MSteps
        Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
        Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,1:nElemsX,1:nElemsY)
#endif
      END DO
    END DO
  END IF

  IF (kk .NE. KCorr) THEN
    DO ii = 2,MSteps
      U(1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
      WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY)
#endif

      CALL FVTimeDerivative(tStage(ii))
      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WCt(1:nVar,1:nElemsX,1:nElemsY)
#endif

    END DO
  END IF

  Up=Ua

#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

END DO

U(1:nVar,1:nElemsX,1:nElemsY) = Ua(MSteps,1:nVar,1:nElemsX,1:nElemsY)
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


#ifdef CENTEREDPRIMITIVE
  Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
  Wp_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
#endif

END DO


CALL FVTimeDerivative(tStage(1))

DO ii = 1,MSteps
  FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)


#ifdef CENTEREDPRIMITIVE
  FWp_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WCt(1:nVar,1:nElemsX,1:nElemsY)
#endif

END DO

DO kk = 1,KCorr
  IF (kk == KCorr) THEN
    ii = MSteps
    Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(1,1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(1,1:nVar,1:nElemsX,1:nElemsY)
#endif

    DO jj = 1,MSteps
      Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
      Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,1:nElemsX,1:nElemsY)
#endif
    END DO
  ELSE
    DO ii = 2,MSteps
      Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(1,1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
      Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(1,1:nVar,1:nElemsX,1:nElemsY)
#endif

      DO jj = 1,MSteps
        Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
        Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,1:nElemsX,1:nElemsY)
#endif
      END DO
    END DO
  END IF

  IF (kk .NE. KCorr) THEN
    DO ii = 2,MSteps
      U(1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
      WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY)
#endif

      CALL FVTimeDerivative(tStage(ii))
      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WCt(1:nVar,1:nElemsX,1:nElemsY)
#endif

    END DO
  END IF

  Up=Ua

#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

END DO

U(1:nVar,1:nElemsX,1:nElemsY) = Ua(MSteps,1:nVar,1:nElemsX,1:nElemsY)
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


#ifdef CENTEREDPRIMITIVE
  Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
  Wp_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
#endif

END DO


CALL FVTimeDerivative(tStage(1))

DO ii = 1,MSteps
  FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)


#ifdef CENTEREDPRIMITIVE
  FWp_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WCt(1:nVar,1:nElemsX,1:nElemsY)
#endif

END DO

DO kk = 1,KCorr
  IF (kk == KCorr) THEN
    ii = MSteps
    Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(1,1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(1,1:nVar,1:nElemsX,1:nElemsY)
#endif

    DO jj = 1,MSteps
      Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
      Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,1:nElemsX,1:nElemsY)
#endif
    END DO
  ELSE
    DO ii = 2,MSteps
      Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(1,1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
      Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(1,1:nVar,1:nElemsX,1:nElemsY)
#endif

      DO jj = 1,MSteps
        Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
        Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY) + dt*thetaDeC(ii,jj)*FWp_WC(jj,1:nVar,1:nElemsX,1:nElemsY)
#endif
      END DO
    END DO
  END IF

  IF (kk .NE. KCorr) THEN
    DO ii = 2,MSteps
      U(1:nVar,1:nElemsX,1:nElemsY) = Ua(ii,1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
      WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(ii,1:nVar,1:nElemsX,1:nElemsY)
#endif

      CALL FVTimeDerivative(tStage(ii))
      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,1:nElemsX,1:nElemsY) = WCt(1:nVar,1:nElemsX,1:nElemsY)
#endif

    END DO
  END IF

  Up=Ua

#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

END DO

U(1:nVar,1:nElemsX,1:nElemsY) = Ua(MSteps,1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(MSteps,1:nVar,1:nElemsX,1:nElemsY)
#endif


!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationByDeC5
!===============================================================================!
!
!
!===============================================================================!
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef CENTEREDPRIMITIVE
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


#ifdef CENTEREDPRIMITIVE
!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: K0_WC
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
REAL            :: initial_guess(1:NRows_WC)
#endif
REAL            :: px,py,dx,dy
!-------------------------------------------------------------------------------!
!*Debug
!-------------------------------------------------------------------------------!
REAL                                :: residual, dxvx, dyvy, partx, party
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

#ifdef CENTEREDPRIMITIVE
old_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)=WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
!*=============================================
!*=============================================
!*=============================================
#endif

K0(1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
K0_WC(1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY) 
#endif

tStage = t 
CALL FVTimeDerivative(tStage)

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
#ifdef CENTEREDPRIMITIVE
WC(1,1:nElemsX,1:nElemsY) = K0_WC(1,1:nElemsX,1:nElemsY) + WCt(1,1:nElemsX,1:nElemsY)*dt 
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
!===============================================================================!
!
!
!
!===============================================================================!
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
!===============================================================================!
!
!
!
!===============================================================================!
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
      rhs_WC(indc)=Wtemp(4)
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
      rhs_WC(indc)=Wtemp(4)
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
      rhs_WC(indc)=Wtemp(4)
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
      rhs_WC(indc)=Wtemp(4)
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
      rhs_WC(indc)=Wtemp(4)
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
      rhs_WC(indc)=Wtemp(4)
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
      rhs_WC(indc)=Wtemp(4)
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
      rhs_WC(indc)=Wtemp(4)
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
#ifdef CENTEREDPRIMITIVE
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


#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC 
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC
#endif

!*NB: BCs are crucial here
USE MOD_Equation           ,ONLY: BoundaryConditions


#ifdef CENTEREDPRIMITIVE
USE MOD_Equation           ,ONLY: Impose_BC_on_WCt
USE MOD_Equation           ,ONLY: Impose_BC_on_WC
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix
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
#ifdef CENTEREDPRIMITIVE
REAL                 :: initial_guess_WC(1:NRows_WC)
#endif

REAL                 :: px,py,dx,dy,ux,vy,part_x,part_y
REAL                 :: FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)
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
#ifdef CENTEREDPRIMITIVE
CALL Put_Matrix_In_Vector(WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_WC(1:NRows_WC))
#endif

DO ii = 1,MSteps
  Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
  Up(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)

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


#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

  !*---------------------------------------------
  !*Impose BCs on p
  !*---------------------------------------------
  IF (kk .NE. 1) THEN
    DO ii=2,MSteps

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

#ifdef CENTEREDPRIMITIVE
    WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
    !*==============================
    CALL FVTimeDerivative(tStage(1)) 

    !*------------------
    !*Add extra splitted contribution
    !*------------------
    !*NB: I have already imposed BCs on Wp
    !*------------------



#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY  
      DO indi=1,nElemsX
        px=First_Derivative_Central_Order2(Wp_WC(1,4,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_WC(1,4,indi,indj-1:indj+1),dy)
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



#ifdef CENTEREDPRIMITIVE
    CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

    DO ii = 1,MSteps
      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)


#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  ELSE !*Other iterations, fluxes have to be computed in all subtimenodes but the first one
    DO ii = 2,MSteps !*NB: The first one does not change
      U(1:nVar,1:nElemsX,1:nElemsY) = Up(ii,1:nVar,1:nElemsX,1:nElemsY)


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

      !*------------------
      !*Add extra splitted contribution
      !*------------------
      !*NB: I have already imposed BCs on Wp
      !*------------------

#ifdef CENTEREDPRIMITIVE
      DO indj=1,nElemsY  
        DO indi=1,nElemsX
          px=First_Derivative_Central_Order2(Wp_WC(ii,4,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_WC(ii,4,indi,indj-1:indj+1),dy)
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




#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)


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


#ifdef CENTEREDPRIMITIVE
    FLUXES_WC = 0.0
#endif

    DO jj = 1,MSteps

      FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) = FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) + thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)


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

#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_WC(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    !*->Linear system for p


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
    CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),Wa_WC(ii,4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
#endif

    !*->Modified update for v
    dx=MESH_DX(1)
    dy=MESH_DX(2)
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------

#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY
      DO indi=1,nElemsX
        Wa_WC(ii,2:3,indi,indj) = Wp_WC(1,2:3,indi,indj) + dt*FLUXES_WC(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_WC(ii,4,indi-1:indi+1,indj)-Wp_WC(ii,4,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_WC(ii,4,indi,indj-1:indj+1)-Wp_WC(ii,4,indi,indj-1:indj+1),dy)
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


#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC 
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC
#endif

!*NB: BCs are crucial here
USE MOD_Equation           ,ONLY: BoundaryConditions


#ifdef CENTEREDPRIMITIVE
USE MOD_Equation           ,ONLY: Impose_BC_on_WCt
USE MOD_Equation           ,ONLY: Impose_BC_on_WC
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix
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
#ifdef CENTEREDPRIMITIVE
REAL                 :: initial_guess_WC(1:NRows_WC)
#endif

REAL                 :: px,py,dx,dy,ux,vy,part_x,part_y
REAL                 :: FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)
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
#ifdef CENTEREDPRIMITIVE
CALL Put_Matrix_In_Vector(WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_WC(1:NRows_WC))
#endif

DO ii = 1,MSteps
  Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
  Up(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)

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


#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

  !*---------------------------------------------
  !*Impose BCs on p
  !*---------------------------------------------
  IF (kk .NE. 1) THEN
    DO ii=2,MSteps

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

#ifdef CENTEREDPRIMITIVE
    WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
    !*==============================
    CALL FVTimeDerivative(tStage(1)) 

    !*------------------
    !*Add extra splitted contribution
    !*------------------
    !*NB: I have already imposed BCs on Wp
    !*------------------



#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY  
      DO indi=1,nElemsX
        px=First_Derivative_Central_Order2(Wp_WC(1,4,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_WC(1,4,indi,indj-1:indj+1),dy)
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



#ifdef CENTEREDPRIMITIVE
    CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

    DO ii = 1,MSteps
      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)


#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  ELSE !*Other iterations, fluxes have to be computed in all subtimenodes but the first one
    DO ii = 2,MSteps !*NB: The first one does not change
      U(1:nVar,1:nElemsX,1:nElemsY) = Up(ii,1:nVar,1:nElemsX,1:nElemsY)


#ifdef CENTEREDPRIMITIVE
      WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

#if(1==1)
      !*Stage dependent max_ro and min_p
      CALL Compute_min_p_max_ro()
#endif

      CALL FVTimeDerivative(tStage(ii))

      !*------------------
      !*Add extra splitted contribution
      !*------------------
      !*NB: I have already imposed BCs on Wp
      !*------------------

#ifdef CENTEREDPRIMITIVE
      DO indj=1,nElemsY  
        DO indi=1,nElemsX
          px=First_Derivative_Central_Order2(Wp_WC(ii,4,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_WC(ii,4,indi,indj-1:indj+1),dy)
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




#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)


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


#ifdef CENTEREDPRIMITIVE
    FLUXES_WC = 0.0
#endif

    DO jj = 1,MSteps

      FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) = FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) + thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)


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

#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_WC(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    !*->Linear system for p


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
    CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),Wa_WC(ii,4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
#endif

    !*->Modified update for v
    dx=MESH_DX(1)
    dy=MESH_DX(2)
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------

#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY
      DO indi=1,nElemsX
        Wa_WC(ii,2:3,indi,indj) = Wp_WC(1,2:3,indi,indj) + dt*FLUXES_WC(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_WC(ii,4,indi-1:indi+1,indj)-Wp_WC(ii,4,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_WC(ii,4,indi,indj-1:indj+1)-Wp_WC(ii,4,indi,indj-1:indj+1),dy)
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
#ifdef CENTEREDPRIMITIVE
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


#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC 
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC
#endif

!*NB: BCs are crucial here
USE MOD_Equation           ,ONLY: BoundaryConditions


#ifdef CENTEREDPRIMITIVE
USE MOD_Equation           ,ONLY: Impose_BC_on_WCt
USE MOD_Equation           ,ONLY: Impose_BC_on_WC
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix
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
#ifdef CENTEREDPRIMITIVE
REAL                 :: initial_guess_WC(1:NRows_WC)
#endif

REAL                 :: px,py,dx,dy,ux,vy,part_x,part_y
REAL                 :: FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)
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
#ifdef CENTEREDPRIMITIVE
CALL Put_Matrix_In_Vector(WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_WC(1:NRows_WC))
#endif

DO ii = 1,MSteps
  Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
  Up(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)

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


#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

  !*---------------------------------------------
  !*Impose BCs on p
  !*---------------------------------------------
  IF (kk .NE. 1) THEN
    DO ii=2,MSteps

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

#ifdef CENTEREDPRIMITIVE
    WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
    !*==============================
    CALL FVTimeDerivative(tStage(1)) 

    !*------------------
    !*Add extra splitted contribution
    !*------------------
    !*NB: I have already imposed BCs on Wp
    !*------------------



#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY  
      DO indi=1,nElemsX
        px=First_Derivative_Central_Order2(Wp_WC(1,4,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_WC(1,4,indi,indj-1:indj+1),dy)
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



#ifdef CENTEREDPRIMITIVE
    CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

    DO ii = 1,MSteps
      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)


#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  ELSE !*Other iterations, fluxes have to be computed in all subtimenodes but the first one
    DO ii = 2,MSteps !*NB: The first one does not change
      U(1:nVar,1:nElemsX,1:nElemsY) = Up(ii,1:nVar,1:nElemsX,1:nElemsY)


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

      !*------------------
      !*Add extra splitted contribution
      !*------------------
      !*NB: I have already imposed BCs on Wp
      !*------------------

#ifdef CENTEREDPRIMITIVE
      DO indj=1,nElemsY  
        DO indi=1,nElemsX
          px=First_Derivative_Central_Order2(Wp_WC(ii,4,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_WC(ii,4,indi,indj-1:indj+1),dy)
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




#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)


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


#ifdef CENTEREDPRIMITIVE
    FLUXES_WC = 0.0
#endif

    DO jj = 1,MSteps

      FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) = FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) + thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)


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

#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_WC(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    !*->Linear system for p


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
    CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),Wa_WC(ii,4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
#endif

    !*->Modified update for v
    dx=MESH_DX(1)
    dy=MESH_DX(2)
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------

#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY
      DO indi=1,nElemsX
        Wa_WC(ii,2:3,indi,indj) = Wp_WC(1,2:3,indi,indj) + dt*FLUXES_WC(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_WC(ii,4,indi-1:indi+1,indj)-Wp_WC(ii,4,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_WC(ii,4,indi,indj-1:indj+1)-Wp_WC(ii,4,indi,indj-1:indj+1),dy)
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
#ifdef CENTEREDPRIMITIVE
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


#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC 
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC
#endif

!*NB: BCs are crucial here
USE MOD_Equation           ,ONLY: BoundaryConditions


#ifdef CENTEREDPRIMITIVE
USE MOD_Equation           ,ONLY: Impose_BC_on_WCt
USE MOD_Equation           ,ONLY: Impose_BC_on_WC
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix
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
#ifdef CENTEREDPRIMITIVE
REAL                 :: initial_guess_WC(1:NRows_WC)
#endif

REAL                 :: px,py,dx,dy,ux,vy,part_x,part_y
REAL                 :: FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)
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
#ifdef CENTEREDPRIMITIVE
CALL Put_Matrix_In_Vector(WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_WC(1:NRows_WC))
#endif

DO ii = 1,MSteps
  Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
  Up(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)

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


#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

  !*---------------------------------------------
  !*Impose BCs on p
  !*---------------------------------------------
  IF (kk .NE. 1) THEN
    DO ii=2,MSteps

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

#ifdef CENTEREDPRIMITIVE
    WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
    !*==============================
    CALL FVTimeDerivative(tStage(1)) 

    !*------------------
    !*Add extra splitted contribution
    !*------------------
    !*NB: I have already imposed BCs on Wp
    !*------------------



#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY  
      DO indi=1,nElemsX
        px=First_Derivative_Central_Order2(Wp_WC(1,4,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_WC(1,4,indi,indj-1:indj+1),dy)
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



#ifdef CENTEREDPRIMITIVE
    CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

    DO ii = 1,MSteps
      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)


#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  ELSE !*Other iterations, fluxes have to be computed in all subtimenodes but the first one
    DO ii = 2,MSteps !*NB: The first one does not change
      U(1:nVar,1:nElemsX,1:nElemsY) = Up(ii,1:nVar,1:nElemsX,1:nElemsY)


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

      !*------------------
      !*Add extra splitted contribution
      !*------------------
      !*NB: I have already imposed BCs on Wp
      !*------------------

#ifdef CENTEREDPRIMITIVE
      DO indj=1,nElemsY  
        DO indi=1,nElemsX
          px=First_Derivative_Central_Order2(Wp_WC(ii,4,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_WC(ii,4,indi,indj-1:indj+1),dy)
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




#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)


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


#ifdef CENTEREDPRIMITIVE
    FLUXES_WC = 0.0
#endif

    DO jj = 1,MSteps

      FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) = FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) + thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)


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

#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_WC(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    !*->Linear system for p


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
    CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),Wa_WC(ii,4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
#endif

    !*->Modified update for v
    dx=MESH_DX(1)
    dy=MESH_DX(2)
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------

#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY
      DO indi=1,nElemsX
        Wa_WC(ii,2:3,indi,indj) = Wp_WC(1,2:3,indi,indj) + dt*FLUXES_WC(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_WC(ii,4,indi-1:indi+1,indj)-Wp_WC(ii,4,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_WC(ii,4,indi,indj-1:indj+1)-Wp_WC(ii,4,indi,indj-1:indj+1),dy)
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
#ifdef CENTEREDPRIMITIVE
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


#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC 
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC
#endif

!*NB: BCs are crucial here
USE MOD_Equation           ,ONLY: BoundaryConditions


#ifdef CENTEREDPRIMITIVE
USE MOD_Equation           ,ONLY: Impose_BC_on_WCt
USE MOD_Equation           ,ONLY: Impose_BC_on_WC
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix
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
#ifdef CENTEREDPRIMITIVE
REAL                 :: initial_guess_WC(1:NRows_WC)
#endif

REAL                 :: px,py,dx,dy,ux,vy,part_x,part_y
REAL                 :: FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)
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
#ifdef CENTEREDPRIMITIVE
CALL Put_Matrix_In_Vector(WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_WC(1:NRows_WC))
#endif

DO ii = 1,MSteps
  Ua(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
  Up(ii,1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)

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


#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

  !*---------------------------------------------
  !*Impose BCs on p
  !*---------------------------------------------
  IF (kk .NE. 1) THEN
    DO ii=2,MSteps

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

#ifdef CENTEREDPRIMITIVE
    WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
    !*==============================
    CALL FVTimeDerivative(tStage(1)) 

    !*------------------
    !*Add extra splitted contribution
    !*------------------
    !*NB: I have already imposed BCs on Wp
    !*------------------



#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY  
      DO indi=1,nElemsX
        px=First_Derivative_Central_Order2(Wp_WC(1,4,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_WC(1,4,indi,indj-1:indj+1),dy)
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



#ifdef CENTEREDPRIMITIVE
    CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

    DO ii = 1,MSteps
      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)


#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  ELSE !*Other iterations, fluxes have to be computed in all subtimenodes but the first one
    DO ii = 2,MSteps !*NB: The first one does not change
      U(1:nVar,1:nElemsX,1:nElemsY) = Up(ii,1:nVar,1:nElemsX,1:nElemsY)


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

      !*------------------
      !*Add extra splitted contribution
      !*------------------
      !*NB: I have already imposed BCs on Wp
      !*------------------

#ifdef CENTEREDPRIMITIVE
      DO indj=1,nElemsY  
        DO indi=1,nElemsX
          px=First_Derivative_Central_Order2(Wp_WC(ii,4,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_WC(ii,4,indi,indj-1:indj+1),dy)
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




#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

      FUp(ii,1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)


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


#ifdef CENTEREDPRIMITIVE
    FLUXES_WC = 0.0
#endif

    DO jj = 1,MSteps

      FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) = FLUXES_U(1:nVar,1:nElemsX,1:nElemsY) + thetaDeC(ii,jj)*FUp(jj,1:nVar,1:nElemsX,1:nElemsY)


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

#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_WC(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    !*->Linear system for p


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
    CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),Wa_WC(ii,4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
#endif

    !*->Modified update for v
    dx=MESH_DX(1)
    dy=MESH_DX(2)
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------

#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY
      DO indi=1,nElemsX
        Wa_WC(ii,2:3,indi,indj) = Wp_WC(1,2:3,indi,indj) + dt*FLUXES_WC(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_WC(ii,4,indi-1:indi+1,indj)-Wp_WC(ii,4,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_WC(ii,4,indi,indj-1:indj+1)-Wp_WC(ii,4,indi,indj-1:indj+1),dy)
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
#ifdef CENTEREDPRIMITIVE
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


#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: K0_WC
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
!-------------------------------------------------------------------------------!
!*Debug
!-------------------------------------------------------------------------------!
REAL                                :: residual, dxvx, dyvy, partx, party
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
#ifdef CENTEREDPRIMITIVE
old_WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)=WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
!*=============================================
!*=============================================
!*=============================================
#endif

K0(1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)
#ifdef CENTEREDPRIMITIVE
K0_WC(1:nVar,1:nElemsX,1:nElemsY) = WC(1:nVar,1:nElemsX,1:nElemsY)
#endif

tStage = t 
CALL FVTimeDerivativePrimitiveOnly(tStage)


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
#ifdef CENTEREDPRIMITIVE
WC(1,1:nElemsX,1:nElemsY) = K0_WC(1,1:nElemsX,1:nElemsY) + WCt(1,1:nElemsX,1:nElemsY)*dt 
#endif


!*=============================

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
#ifdef CENTEREDPRIMITIVE
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


#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC !W_WC^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC  !F(W_WC^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC  !W_WC^(k)
#endif

!*NB: BCs are crucial here
USE MOD_Equation           ,ONLY: BoundaryConditions

#ifdef CENTEREDPRIMITIVE
USE MOD_Equation           ,ONLY: Impose_BC_on_WCt
USE MOD_Equation           ,ONLY: Impose_BC_on_WC
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix
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
#ifdef CENTEREDPRIMITIVE
REAL                 :: initial_guess_WC(1:NRows_WC)
#endif
REAL                 :: px,py,dx,dy,ux,vy,part_x,part_y
REAL                 :: FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)
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

#ifdef CENTEREDPRIMITIVE
CALL Put_Matrix_In_Vector(WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_WC(1:NRows_WC))
#endif

DO ii = 1,MSteps

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
#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

  !*---------------------------------------------
  !*Impose BCs on p
  !*---------------------------------------------
  IF (kk .NE. 1) THEN
    DO ii=2,MSteps
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
#ifdef CENTEREDPRIMITIVE
    WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
    !*==============================
    CALL FVTimeDerivativePrimitiveOnly(tStage(1)) 

    !*------------------
    !*Add extra splitted contribution
    !*------------------
    !*NB: I have already imposed BCs on Wp
    !*------------------


#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY  
      DO indi=1,nElemsX
        px=First_Derivative_Central_Order2(Wp_WC(1,4,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_WC(1,4,indi,indj-1:indj+1),dy)
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



#ifdef CENTEREDPRIMITIVE
    CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

    DO ii = 1,MSteps

#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  ELSE !*Other iterations, fluxes have to be computed in all subtimenodes but the first one
    DO ii = 2,MSteps !*NB: The first one does not change

#ifdef CENTEREDPRIMITIVE
      WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

      !*======================================
      !*NB: HERE I NEED TO CALL FV TIME DERIVATIVE OF THE PRIMITIVE VARIABLES ONLY
      !*======================================
      CALL FVTimeDerivativePrimitiveOnly(tStage(ii)) 

      !*------------------
      !*Add extra splitted contribution
      !*------------------
      !*NB: I have already imposed BCs on Wp
      !*------------------

#ifdef CENTEREDPRIMITIVE
      DO indj=1,nElemsY  
        DO indi=1,nElemsX
          px=First_Derivative_Central_Order2(Wp_WC(ii,4,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_WC(ii,4,indi,indj-1:indj+1),dy)
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



#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
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
#ifdef CENTEREDPRIMITIVE
    FLUXES_WC = 0.0
#endif

    DO jj = 1,MSteps

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
#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_WC(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    !*->Linear system for p


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
    CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),Wa_WC(ii,4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
#endif

    !*->Modified update for v
    dx=MESH_DX(1)
    dy=MESH_DX(2)
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------

#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY
      DO indi=1,nElemsX
        Wa_WC(ii,2:3,indi,indj) = Wp_WC(1,2:3,indi,indj) + dt*FLUXES_WC(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_WC(ii,4,indi-1:indi+1,indj)-Wp_WC(ii,4,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_WC(ii,4,indi,indj-1:indj+1)-Wp_WC(ii,4,indi,indj-1:indj+1),dy)
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

#ifdef CENTEREDPRIMITIVE
  CALL Impose_BC_on_WC(Wa_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),tStage(ii))
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
#ifdef CENTEREDPRIMITIVE
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


#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC !W_WC^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC  !F(W_WC^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC  !W_WC^(k)
#endif

!*NB: BCs are crucial here
USE MOD_Equation           ,ONLY: BoundaryConditions

#ifdef CENTEREDPRIMITIVE
USE MOD_Equation           ,ONLY: Impose_BC_on_WCt
USE MOD_Equation           ,ONLY: Impose_BC_on_WC
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix
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
#ifdef CENTEREDPRIMITIVE
REAL                 :: initial_guess_WC(1:NRows_WC)
#endif
REAL                 :: px,py,dx,dy,ux,vy,part_x,part_y
REAL                 :: FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)
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

#ifdef CENTEREDPRIMITIVE
CALL Put_Matrix_In_Vector(WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_WC(1:NRows_WC))
#endif

DO ii = 1,MSteps

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
#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

  !*---------------------------------------------
  !*Impose BCs on p
  !*---------------------------------------------
  IF (kk .NE. 1) THEN
    DO ii=2,MSteps
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
#ifdef CENTEREDPRIMITIVE
    WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
    !*==============================
    CALL FVTimeDerivativePrimitiveOnly(tStage(1)) 

    !*------------------
    !*Add extra splitted contribution
    !*------------------
    !*NB: I have already imposed BCs on Wp
    !*------------------


#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY  
      DO indi=1,nElemsX
        px=First_Derivative_Central_Order2(Wp_WC(1,4,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_WC(1,4,indi,indj-1:indj+1),dy)
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



#ifdef CENTEREDPRIMITIVE
    CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

    DO ii = 1,MSteps

#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  ELSE !*Other iterations, fluxes have to be computed in all subtimenodes but the first one
    DO ii = 2,MSteps !*NB: The first one does not change

#ifdef CENTEREDPRIMITIVE
      WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

      !*======================================
      !*NB: HERE I NEED TO CALL FV TIME DERIVATIVE OF THE PRIMITIVE VARIABLES ONLY
      !*======================================
      CALL FVTimeDerivativePrimitiveOnly(tStage(ii)) 

      !*------------------
      !*Add extra splitted contribution
      !*------------------
      !*NB: I have already imposed BCs on Wp
      !*------------------

#ifdef CENTEREDPRIMITIVE
      DO indj=1,nElemsY  
        DO indi=1,nElemsX
          px=First_Derivative_Central_Order2(Wp_WC(ii,4,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_WC(ii,4,indi,indj-1:indj+1),dy)
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



#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
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
#ifdef CENTEREDPRIMITIVE
    FLUXES_WC = 0.0
#endif

    DO jj = 1,MSteps

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
#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_WC(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    !*->Linear system for p


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
    CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),Wa_WC(ii,4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
#endif

    !*->Modified update for v
    dx=MESH_DX(1)
    dy=MESH_DX(2)
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------

#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY
      DO indi=1,nElemsX
        Wa_WC(ii,2:3,indi,indj) = Wp_WC(1,2:3,indi,indj) + dt*FLUXES_WC(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_WC(ii,4,indi-1:indi+1,indj)-Wp_WC(ii,4,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_WC(ii,4,indi,indj-1:indj+1)-Wp_WC(ii,4,indi,indj-1:indj+1),dy)
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

#ifdef CENTEREDPRIMITIVE
  CALL Impose_BC_on_WC(Wa_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),tStage(ii))
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
#ifdef CENTEREDPRIMITIVE
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


#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC !W_WC^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC  !F(W_WC^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC  !W_WC^(k)
#endif

!*NB: BCs are crucial here
USE MOD_Equation           ,ONLY: BoundaryConditions

#ifdef CENTEREDPRIMITIVE
USE MOD_Equation           ,ONLY: Impose_BC_on_WCt
USE MOD_Equation           ,ONLY: Impose_BC_on_WC
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix
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
#ifdef CENTEREDPRIMITIVE
REAL                 :: initial_guess_WC(1:NRows_WC)
#endif
REAL                 :: px,py,dx,dy,ux,vy,part_x,part_y
REAL                 :: FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)
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

#ifdef CENTEREDPRIMITIVE
CALL Put_Matrix_In_Vector(WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_WC(1:NRows_WC))
#endif

DO ii = 1,MSteps

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
#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

  !*---------------------------------------------
  !*Impose BCs on p
  !*---------------------------------------------
  IF (kk .NE. 1) THEN
    DO ii=2,MSteps
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
#ifdef CENTEREDPRIMITIVE
    WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
    !*==============================
    CALL FVTimeDerivativePrimitiveOnly(tStage(1)) 

    !*------------------
    !*Add extra splitted contribution
    !*------------------
    !*NB: I have already imposed BCs on Wp
    !*------------------


#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY  
      DO indi=1,nElemsX
        px=First_Derivative_Central_Order2(Wp_WC(1,4,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_WC(1,4,indi,indj-1:indj+1),dy)
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



#ifdef CENTEREDPRIMITIVE
    CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

    DO ii = 1,MSteps

#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  ELSE !*Other iterations, fluxes have to be computed in all subtimenodes but the first one
    DO ii = 2,MSteps !*NB: The first one does not change

#ifdef CENTEREDPRIMITIVE
      WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

      !*======================================
      !*NB: HERE I NEED TO CALL FV TIME DERIVATIVE OF THE PRIMITIVE VARIABLES ONLY
      !*======================================
      CALL FVTimeDerivativePrimitiveOnly(tStage(ii)) 

      !*------------------
      !*Add extra splitted contribution
      !*------------------
      !*NB: I have already imposed BCs on Wp
      !*------------------

#ifdef CENTEREDPRIMITIVE
      DO indj=1,nElemsY  
        DO indi=1,nElemsX
          px=First_Derivative_Central_Order2(Wp_WC(ii,4,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_WC(ii,4,indi,indj-1:indj+1),dy)
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



#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
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
#ifdef CENTEREDPRIMITIVE
    FLUXES_WC = 0.0
#endif

    DO jj = 1,MSteps

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
#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_WC(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    !*->Linear system for p


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
    CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),Wa_WC(ii,4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
#endif

    !*->Modified update for v
    dx=MESH_DX(1)
    dy=MESH_DX(2)
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------

#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY
      DO indi=1,nElemsX
        Wa_WC(ii,2:3,indi,indj) = Wp_WC(1,2:3,indi,indj) + dt*FLUXES_WC(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_WC(ii,4,indi-1:indi+1,indj)-Wp_WC(ii,4,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_WC(ii,4,indi,indj-1:indj+1)-Wp_WC(ii,4,indi,indj-1:indj+1),dy)
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

#ifdef CENTEREDPRIMITIVE
  CALL Impose_BC_on_WC(Wa_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),tStage(ii))
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
#ifdef CENTEREDPRIMITIVE
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


#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC !W_WC^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC  !F(W_WC^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC  !W_WC^(k)
#endif

!*NB: BCs are crucial here
USE MOD_Equation           ,ONLY: BoundaryConditions

#ifdef CENTEREDPRIMITIVE
USE MOD_Equation           ,ONLY: Impose_BC_on_WCt
USE MOD_Equation           ,ONLY: Impose_BC_on_WC
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix
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
#ifdef CENTEREDPRIMITIVE
REAL                 :: initial_guess_WC(1:NRows_WC)
#endif
REAL                 :: px,py,dx,dy,ux,vy,part_x,part_y
REAL                 :: FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)
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

#ifdef CENTEREDPRIMITIVE
CALL Put_Matrix_In_Vector(WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_WC(1:NRows_WC))
#endif

DO ii = 1,MSteps

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
#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

  !*---------------------------------------------
  !*Impose BCs on p
  !*---------------------------------------------
  IF (kk .NE. 1) THEN
    DO ii=2,MSteps
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
#ifdef CENTEREDPRIMITIVE
    WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
    !*==============================
    CALL FVTimeDerivativePrimitiveOnly(tStage(1)) 

    !*------------------
    !*Add extra splitted contribution
    !*------------------
    !*NB: I have already imposed BCs on Wp
    !*------------------


#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY  
      DO indi=1,nElemsX
        px=First_Derivative_Central_Order2(Wp_WC(1,4,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_WC(1,4,indi,indj-1:indj+1),dy)
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



#ifdef CENTEREDPRIMITIVE
    CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

    DO ii = 1,MSteps

#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  ELSE !*Other iterations, fluxes have to be computed in all subtimenodes but the first one
    DO ii = 2,MSteps !*NB: The first one does not change

#ifdef CENTEREDPRIMITIVE
      WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

      !*======================================
      !*NB: HERE I NEED TO CALL FV TIME DERIVATIVE OF THE PRIMITIVE VARIABLES ONLY
      !*======================================
      CALL FVTimeDerivativePrimitiveOnly(tStage(ii)) 

      !*------------------
      !*Add extra splitted contribution
      !*------------------
      !*NB: I have already imposed BCs on Wp
      !*------------------

#ifdef CENTEREDPRIMITIVE
      DO indj=1,nElemsY  
        DO indi=1,nElemsX
          px=First_Derivative_Central_Order2(Wp_WC(ii,4,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_WC(ii,4,indi,indj-1:indj+1),dy)
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



#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
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
#ifdef CENTEREDPRIMITIVE
    FLUXES_WC = 0.0
#endif

    DO jj = 1,MSteps

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
#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_WC(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    !*->Linear system for p


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
    CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),Wa_WC(ii,4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
#endif

    !*->Modified update for v
    dx=MESH_DX(1)
    dy=MESH_DX(2)
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------

#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY
      DO indi=1,nElemsX
        Wa_WC(ii,2:3,indi,indj) = Wp_WC(1,2:3,indi,indj) + dt*FLUXES_WC(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_WC(ii,4,indi-1:indi+1,indj)-Wp_WC(ii,4,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_WC(ii,4,indi,indj-1:indj+1)-Wp_WC(ii,4,indi,indj-1:indj+1),dy)
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

#ifdef CENTEREDPRIMITIVE
  CALL Impose_BC_on_WC(Wa_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),tStage(ii))
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
#ifdef CENTEREDPRIMITIVE
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


#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC !W_WC^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC  !F(W_WC^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC  !W_WC^(k)
#endif

!*NB: BCs are crucial here
USE MOD_Equation           ,ONLY: BoundaryConditions


#ifdef CENTEREDPRIMITIVE
USE MOD_Equation           ,ONLY: Impose_BC_on_WCt
USE MOD_Equation           ,ONLY: Impose_BC_on_WC
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix
#endif

USE MOD_JacobiIteration,    ONLY: jacobi_2
USE MOD_IterativeLinearSolver,ONLY: pmgmres_ilu_cr


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
#ifdef CENTEREDPRIMITIVE
REAL                 :: initial_guess_WC(1:NRows_WC)
#endif
REAL                 :: px,py,dx,dy,ux,vy,part_x,part_y
REAL                 :: FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)
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

#ifdef CENTEREDPRIMITIVE
CALL Put_Matrix_In_Vector(WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_WC(1:NRows_WC))
#endif

DO ii = 1,MSteps

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
#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

  !*---------------------------------------------
  !*Impose BCs on p
  !*---------------------------------------------
  IF (kk .NE. 1) THEN
    DO ii=2,MSteps
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
#ifdef CENTEREDPRIMITIVE
    WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
    !*==============================
    CALL FVTimeDerivativePrimitiveOnly(tStage(1)) 

    !*------------------
    !*Add extra splitted contribution
    !*------------------
    !*NB: I have already imposed BCs on Wp
    !*------------------


#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY  
      DO indi=1,nElemsX
        px=First_Derivative_Central_Order2(Wp_WC(1,4,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_WC(1,4,indi,indj-1:indj+1),dy)
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



#ifdef CENTEREDPRIMITIVE
    CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

    DO ii = 1,MSteps

#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  ELSE !*Other iterations, fluxes have to be computed in all subtimenodes but the first one
    DO ii = 2,MSteps !*NB: The first one does not change


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

      !*------------------
      !*Add extra splitted contribution
      !*------------------
      !*NB: I have already imposed BCs on Wp
      !*------------------
#ifdef CENTEREDPRIMITIVE
      DO indj=1,nElemsY  
        DO indi=1,nElemsX
          px=First_Derivative_Central_Order2(Wp_WC(ii,4,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_WC(ii,4,indi,indj-1:indj+1),dy)
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



#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
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
#ifdef CENTEREDPRIMITIVE
    FLUXES_WC = 0.0
#endif

    DO jj = 1,MSteps

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
#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_WC(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    !*->Linear system for p

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
    CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),Wa_WC(ii,4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))

#endif

    !*->Modified update for v
    dx=MESH_DX(1)
    dy=MESH_DX(2)


#ifdef CENTEREDPRIMITIVE
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
    DO indj=1,nElemsY
      DO indi=1,nElemsX
        Wa_WC(ii,2:3,indi,indj) = Wp_WC(1,2:3,indi,indj) + dt*FLUXES_WC(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_WC(ii,4,indi-1:indi+1,indj)-Wp_WC(ii,4,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_WC(ii,4,indi,indj-1:indj+1)-Wp_WC(ii,4,indi,indj-1:indj+1),dy)
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
#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(1,1:nVar,1:nElemsX,1:nElemsY)
#endif

CALL FVTimeDerivativeConservedOnly(t)

K0(1:nVar,1:nElemsX,1:nElemsY)=Ut(1:nVar,1:nElemsX,1:nElemsY)


!*====================
!*tn+1
!*====================
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
#ifdef CENTEREDPRIMITIVE
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


#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC !W_WC^(k-1)
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC  !F(W_WC^(k-1))
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC  !W_WC^(k)
#endif

!*NB: BCs are crucial here
USE MOD_Equation           ,ONLY: BoundaryConditions


#ifdef CENTEREDPRIMITIVE
USE MOD_Equation           ,ONLY: Impose_BC_on_WCt
USE MOD_Equation           ,ONLY: Impose_BC_on_WC
USE MOD_FiniteVolume2D     ,ONLY: Put_Matrix_In_Vector
USE MOD_FiniteVolume2D     ,ONLY: Put_Vector_In_Matrix
#endif

USE MOD_JacobiIteration,    ONLY: jacobi_2
USE MOD_IterativeLinearSolver,ONLY: pmgmres_ilu_cr


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
#ifdef CENTEREDPRIMITIVE
REAL                 :: initial_guess_WC(1:NRows_WC)
#endif
REAL                 :: px,py,dx,dy,ux,vy,part_x,part_y
REAL                 :: FLUXES_U(1:nVar,1:nElemsX,1:nElemsY)
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

#ifdef CENTEREDPRIMITIVE
CALL Put_Matrix_In_Vector(WC(4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),initial_guess_WC(1:NRows_WC))
#endif

DO ii = 1,MSteps

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
#ifdef CENTEREDPRIMITIVE
  Wp_WC=Wa_WC
#endif

  !*---------------------------------------------
  !*Impose BCs on p
  !*---------------------------------------------
  IF (kk .NE. 1) THEN
    DO ii=2,MSteps
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
#ifdef CENTEREDPRIMITIVE
    WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif
    !*==============================
    CALL FVTimeDerivativePrimitiveOnly(tStage(1)) 

    !*------------------
    !*Add extra splitted contribution
    !*------------------
    !*NB: I have already imposed BCs on Wp
    !*------------------


#ifdef CENTEREDPRIMITIVE
    DO indj=1,nElemsY  
      DO indi=1,nElemsX
        px=First_Derivative_Central_Order2(Wp_WC(1,4,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wp_WC(1,4,indi,indj-1:indj+1),dy)
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



#ifdef CENTEREDPRIMITIVE
    CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
#endif

    DO ii = 1,MSteps

#ifdef CENTEREDPRIMITIVE
      FWp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    END DO
  ELSE !*Other iterations, fluxes have to be computed in all subtimenodes but the first one
    DO ii = 2,MSteps !*NB: The first one does not change


#ifdef CENTEREDPRIMITIVE
      WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(ii,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

#if(1==1)
      !*STAGE DEPENDENT max_ro and min_p ACTIVE
      CALL Compute_min_p_max_ro()
#endif

      CALL FVTimeDerivativePrimitiveOnly(tStage(ii))

      !*------------------
      !*Add extra splitted contribution
      !*------------------
      !*NB: I have already imposed BCs on Wp
      !*------------------
#ifdef CENTEREDPRIMITIVE
      DO indj=1,nElemsY  
        DO indi=1,nElemsX
          px=First_Derivative_Central_Order2(Wp_WC(ii,4,indi-1:indi+1,indj),dx)
          py=First_Derivative_Central_Order2(Wp_WC(ii,4,indi,indj-1:indj+1),dy)
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



#ifdef CENTEREDPRIMITIVE
      CALL Impose_BC_on_WCt() !*NB: It is important to impose BCs on the fluxes because they need to be differentiated
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
#ifdef CENTEREDPRIMITIVE
    FLUXES_WC = 0.0
#endif

    DO jj = 1,MSteps

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
#ifdef CENTEREDPRIMITIVE
    Wa_WC(ii,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) = Wp_WC(1,1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) + dt*FLUXES_WC(1,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
#endif

    !*->Linear system for p

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
    CALL Put_Vector_In_Matrix(sol_WC(1:NRows_WC),Wa_WC(ii,4,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))

#endif

    !*->Modified update for v
    dx=MESH_DX(1)
    dy=MESH_DX(2)


#ifdef CENTEREDPRIMITIVE
    !*----------------------
    !*Update velocity with extra explicit contribution <=============
    !*----------------------
    DO indj=1,nElemsY
      DO indi=1,nElemsX
        Wa_WC(ii,2:3,indi,indj) = Wp_WC(1,2:3,indi,indj) + dt*FLUXES_WC(2:3,indi,indj)
        px=First_Derivative_Central_Order2(Wa_WC(ii,4,indi-1:indi+1,indj)-Wp_WC(ii,4,indi-1:indi+1,indj),dx)
        py=First_Derivative_Central_Order2(Wa_WC(ii,4,indi,indj-1:indj+1)-Wp_WC(ii,4,indi,indj-1:indj+1),dy)
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
#ifdef CENTEREDPRIMITIVE
WC(1:nVar,1:nElemsX,1:nElemsY) = Wa_WC(1,1:nVar,1:nElemsX,1:nElemsY)
#endif

CALL FVTimeDerivativeConservedOnly(t)

K0(1:nVar,1:nElemsX,1:nElemsY)=Ut(1:nVar,1:nElemsX,1:nElemsY)


!*====================
!*tn+1
!*====================
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
!===============================================================================!
!
!
!
!===============================================================================!
!===============================================================================!
!
!
!
!===============================================================================!
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
      rhs_WC(indc)=Wtemp(4)
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
      rhs_WC(indc)=Wtemp(4)
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
      rhs_WC(indc)=Wtemp(4)
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
      rhs_WC(indc)=Wtemp(4)
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
      rhs_WC(indc)=Wtemp(4)
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
      rhs_WC(indc)=Wtemp(4)
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
      rhs_WC(indc)=Wtemp(4)
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
      rhs_WC(indc)=Wtemp(4)
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
#endif
!===============================================================================!
END MODULE MOD_TimeDiscretization
!-------------------------------------------------------------------------------!
