!===============================================================================!
MODULE MOD_Parameters
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE InitializeParameters
  MODULE PROCEDURE InitializeParameters
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: InitializeParameters
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
SUBROUTINE InitializeParameters()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: PI
USE MOD_FiniteVolume2D_vars,ONLY: CFL
USE MOD_FiniteVolume2D_vars,ONLY: TEnd
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructionFix
USE MOD_FiniteVolume2D_vars,ONLY: timescheme
USE MOD_FiniteVolume2D_vars,ONLY: WhichRiemannSolver
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
USE MOD_FiniteVolume2D_vars,ONLY: WhichRiemannSolverPrimitiveSystem
#endif
USE MOD_FiniteVolume2D_vars,ONLY: WhichOutput
USE MOD_FiniteVolume2D_vars,ONLY: nOutputFiles
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsType
USE MOD_FiniteVolume2D_vars,ONLY: VarNameVisu
USE MOD_FiniteVolume2D_vars,ONLY: GravitationalPotentialFlag
USE MOD_FiniteVolume2D_vars,ONLY: maxTimeSteps
USE MOD_FiniteVolume2D_vars,ONLY: MStepsMax
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructedVariable
USE MOD_FiniteVolume2D_vars,ONLY: WhichSpeedEstimateForDtComputation
USE MOD_FiniteVolume2D_vars,ONLY: relaxed_dt
USE MOD_FiniteVolume2D_vars,ONLY: NRelaxedTimesteps
#ifdef SW
USE MOD_FiniteVolume2D_vars,ONLY: Gravity
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
#endif
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
USE MOD_FiniteVolume2D_vars,ONLY: VarNameVisu_W
#endif
#ifdef ACTIVEFLUX
USE MOD_FiniteVolume2D_vars,ONLY: AF_PostProcessing
#endif
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#ifdef IMEX
USE MOD_FiniteVolume2D_vars,ONLY: SAFETY_K_TO_AVOID_HYPERBOLICITY_LOSS
#endif
#endif
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
USE MOD_FiniteVolume2D_vars,ONLY: nVar
#ifdef MULTIFLUID
USE MOD_FiniteVolume2D_vars,ONLY: Gmm1
USE MOD_FiniteVolume2D_vars,ONLY: Gmm2
USE MOD_FiniteVolume2D_vars,ONLY: pinf1
USE MOD_FiniteVolume2D_vars,ONLY: pinf2
USE MOD_FiniteVolume2D_vars,ONLY: N_Troubled_Neighbors
#endif
USE MOD_FiniteVolume2D_vars,ONLY: DiscontinuityLocation
USE MOD_FiniteVolume2D_vars,ONLY: Lambda_Lax
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
#if defined(CENTEREDPRIMITIVE) && defined(PATHCONSERVATIVESHOCKDETECTION)
USE MOD_FiniteVolume2D_vars,ONLY: K_coefficient_PCSD
USE MOD_FiniteVolume2D_vars,ONLY: WhichVariableForLimiting_PCSD
#endif
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
USE MOD_FiniteVolume2D_vars,ONLY: maxLenghtJacobiCounter
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem
#endif
#endif
USE MOD_FiniteVolume2D_vars,ONLY: Ma_shock
USE MOD_FiniteVolume2D_vars,ONLY: Ma_vortex
USE MOD_FiniteVolume2D_vars,ONLY: R_EOS
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
CHARACTER(LEN=255) :: ErrorMessage
INTEGER            :: iarg, nargs, fake_integer
REAL               :: iarg_real
CHARACTER(len=32)  :: arg
CHARACTER(LEN=80)  :: NameTest
INTEGER            :: index
!-------------------------------------------------------------------------------!

PRINT*, "--------------------------"
PRINT*, "Initializing parameters   "
PRINT*, "--------------------------"

InitialCondition = 100  

index=0
nargs = command_argument_COUNT()
IF (nargs > index) THEN
   CALL get_command_ARGUMENT(index+1, arg)
   READ(arg, *) iarg
   InitialCondition = iarg
END IF

SELECT CASE(InitialCondition)
!*------------------------------------------
!*NB: Due to the flag, a test for SW can have the same number as a test for Euler
!*But try to avoid it
!*------------------------------------------
!*Boudary conditions are
!*B, R, U, L
!*1 periodic
!*2 transmissive
!*3 inflow
!*4 outflow
!*5 wall
!*------------------------------------------
#ifdef SW
  !*------------------------------------------
  !*[1] Unsteady smooth vortex for SW
  !*------------------------------------------
  CASE(1) !*UNSTEADY SMOOTH VORTEX
    NameTest="Unsteady smooth vortex for SW (NB: to be run with SW flag)"
    TEnd    = 0.1
    Gravity = 9.81
    Kappa   = 0.5*Gravity
    Gmm     = 2.0
    nElemsX = 512
    nElemsY = nElemsX
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/3.0,3.0/)
    BoundaryConditionsType = (/1,1,1,1/) !*PERIODIC BCs
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#else
  !*------------------------------------------
  !*[2] Steady isentropic smooth vortex
  !*------------------------------------------
  CASE(2) 
    NameTest="Steady isentropic vortex"
    TEnd    = 0.1
    Gmm     = 1.4
    nElemsX = 120
    nElemsY = nElemsX
    MESH_X0 = (/-10.0,-10.0/)
    MESH_X1 = (/10.0,10.0/)
    BoundaryConditionsType = (/1,1,1,1/) !*PERIODIC BCs
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
  !*------------------------------------------
  !*[3] Unsteady isentropic smooth vortex
  !*------------------------------------------
  CASE(3) 
    NameTest="Unsteady isentropic vortex"
    TEnd    = 0.1
    Gmm     = 1.4
    nElemsX = 120
    nElemsY = nElemsX
    MESH_X0 = (/-10.0,-10.0/)
    MESH_X1 = (/10.0,10.0/)
    BoundaryConditionsType = (/1,1,1,1/) !*PERIODIC BCs
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0     
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0  
  !*------------------------------------------
  !*[4] Advection of smooth density
  !*------------------------------------------
  CASE(4) 
    NameTest="Advection of smooth density"
    TEnd    = 0.1
    Gmm     = 1.4
    nElemsX = 120
    nElemsY = nElemsX
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/1.0,1.0/)
    BoundaryConditionsType = (/1,1,1,1/) !*PERIODIC BCs
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
  !*------------------------------------------
  !*[5] 1D advection of smooth density in X direction
  !*------------------------------------------
  CASE(5) 
    NameTest="1D advection of smooth density in X direction"
    TEnd    = 0.1
    Gmm     = 1.4
    nElemsX = 120
    nElemsY = nElemsX
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/1.0,1.0/)
    BoundaryConditionsType = (/1,1,1,1/) !*PERIODIC BCs
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0

  CASE(151) !*Advection of sin^4(pix) wave 
    !*=================================
    !*TVD Fluxes for the High-Order ADER Schemes, Toro, Titarev
    !*=================================

    NameTest="1D advection of sin4 density in X direction"
    TEnd    = 2.0
    Gmm     = 1.4
    nElemsX = 120
    nElemsY = nElemsX
    MESH_X0 = (/-1.0,0.0/)
    MESH_X1 = (/ 1.0,1.0/)
    BoundaryConditionsType = (/1,1,1,1/) !*PERIODIC BCs
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0



  !*------------------------------------------
  !*[6] 1D advection of smooth density in Y direction
  !*------------------------------------------
  CASE(6) 
    NameTest="1D advection of smooth density in Y direction"
    TEnd    = 0.1
    Gmm     = 1.4
    nElemsX = 120
    nElemsY = nElemsX
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/1.0,1.0/)
    BoundaryConditionsType = (/1,1,1,1/) !*PERIODIC BCs
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
  !*------------------------------------------
  !*[7] Sod - Explosion problem
  !*------------------------------------------
  CASE(7) 
    NameTest="Sod - Explosion problem"
    TEnd    = 0.25
    Gmm     = 1.4
    nElemsX = 120
    nElemsY = nElemsX
    MESH_X0 = (/-1.0,-1.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/2,2,2,2/) !*Reflecting
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 0
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif

  !*------------------------------------------
  !*[77] Sod - Explosion problem
  !*------------------------------------------
  CASE(77) 
    NameTest="Sod - Explosion problem"
    TEnd    = 0.25
    Gmm     = 1.4
    nElemsX = 120
    nElemsY = nElemsX
    MESH_X0 = (/-1.0,-1.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/2,2,2,2/) !*Reflecting
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif


  !*------------------------------------------
  !*[-7] Sod - Explosion problem low Mach experiment
  !*------------------------------------------
  CASE(-7) 
    NameTest="Sod - Explosion problem low Mach experiment"
    TEnd    = 0.05
    Gmm     = 1.4
    nElemsX = 120
    nElemsY = nElemsX
    MESH_X0 = (/-1.0,-1.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/2,2,2,2/) !*Reflecting
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 0
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif

  !*------------------------------------------
  !*[-70] Sod - Explosion problem low Mach experiment
  !*------------------------------------------
  CASE(-70) 
    NameTest="Sod - Explosion problem low Mach experiment with time step relaxed"
    TEnd    = 0.05
    Gmm     = 1.4
    nElemsX = 120
    nElemsY = nElemsX
    MESH_X0 = (/-1.0,-1.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/2,2,2,2/) !*Reflecting
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif


  !*------------------------------------------
  !*[8] Sod - 1D in X direction
  !*------------------------------------------
  CASE(8) 
    NameTest="Sod - 1D in X direction"
    TEnd    = 0.2
    Gmm     = 1.4
    nElemsX = 200
    nElemsY = 3
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/2,2,2,2/) !*Reflecting
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 0
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif

  !*------------------------------------------
  !*[88] Sod - 1D in X direction
  !*------------------------------------------
  CASE(88) 
    NameTest="Sod - 1D in X direction"
    TEnd    = 0.2
    Gmm     = 1.4
    nElemsX = 200
    nElemsY = 3
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/2,2,2,2/) !*Reflecting
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif

  !*------------------------------------------
  !*[-8] Sod - 1D in X direction
  !*------------------------------------------
  CASE(-8) 
    NameTest="Sod - 1D in X direction low Mach experiment"
    TEnd    = 0.05
    Gmm     = 1.4
    nElemsX = 200
    nElemsY = 3
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/2,2,2,2/) !*Reflecting
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 0
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif

  !*------------------------------------------
  !*[-80] Sod - 1D in X direction
  !*------------------------------------------
  CASE(-80) 
    NameTest="Sod - 1D in X direction low Mach experiment with time step relaxed"
    TEnd    = 0.05
    Gmm     = 1.4
    nElemsX = 200
    nElemsY = 3
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/2,2,2,2/) !*Reflecting
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif

  !*------------------------------------------
  !*[9] Sod - 1D in Y direction
  !*------------------------------------------
  CASE(9) 
    NameTest="Sod - 1D in Y direction"
    TEnd    = 0.2
    Gmm     = 1.4
    nElemsX = 120
    nElemsY = nElemsX
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/2,2,2,2/) !*Reflecting
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif

  !*------------------------------------------
  !*[10] Implosion problem
  !*------------------------------------------
  CASE(10) 
    NameTest="Implosion problem"
    TEnd    = 2.5 !0.01 !2.5
    Gmm     = 1.4
    nElemsX = 250
    nElemsY = nElemsX
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 0.3, 0.3/)
    BoundaryConditionsType = (/5,5,5,5/) !*Reflecting
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif

  !*------------------------------------------
  !*[11] 2D-Riemann problem
  !*------------------------------------------
  CASE(11) 
    NameTest="2D-Riemann problem"
    TEnd    = 1.0
    Gmm     = 1.4
    nElemsX = 1000
    nElemsY = nElemsX
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 1.2, 1.2/)
    BoundaryConditionsType = (/2,2,2,2/) !*Transmissive (free)
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif


  !*------------------------------------------
  !*[12] Rayleigh-Taylor instability
  !*------------------------------------------
  CASE(12) 
    NameTest="Rayleigh-Taylor instability"
    TEnd    = 2.95
    Gmm     = 5.0/3.0
    nElemsX = 256
    nElemsY = 1024
    MESH_X0 = (/ 0.0 , 0.0/)
    MESH_X1 = (/ 0.25, 1.0/)
    BoundaryConditionsType = (/3,5,3,5/) !*Transmissive (free)
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef MULTIFLUID
    PrimRefState1 = (/2.0, 0.0, 0.0, 1.0, 0.0/)
    PrimRefState3 = (/1.0, 0.0, 0.0, 2.5, 0.0/)
#else
    PrimRefState1 = (/2.0, 0.0, 0.0, 1.0/)
    PrimRefState3 = (/1.0, 0.0, 0.0, 2.5/)
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 1       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif

  !*------------------------------------------
  !*[13] Shock-vortex interaction
  !*------------------------------------------
  CASE(13) 
    NameTest="Shock-vortex interaction"
    TEnd    = 0.69
    Gmm     = 1.4
    nElemsX = 800
    nElemsY = 401
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 2.0, 1.0/)
    Ma_shock  = 1.5
    Ma_vortex = 0.9
    R_EOS     = 287.0
#ifdef MULTIFLUID
    PrimRefState4 = (/1.0, SQRT(Gmm)*Ma_shock,0.0, 1.0,0.0/)
#else
    PrimRefState4 = (/1.0, SQRT(Gmm)*Ma_shock,0.0, 1.0/)
#endif
    BoundaryConditionsType = (/2,2,2,3/) !*Transmissive all but left
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif



  CASE(-300) ! Shock-turbulence interaction Shu-Osher
    !*Efficient implementation of essentially non-oscillatory shock capturing schemes, II, Shu, Osher
    TEnd    = 1.8
    Gmm     = 1.4
#ifdef MULTIFLUID
    PrimRefState4 = (/3.857143, 2.629369,0.0, 10.333333,0.0/)
#else
    PrimRefState4 = (/3.857143, 2.629369,0.0, 10.333333/)
#endif
    nElemsX = 200
    nElemsY = 5
    MESH_X0 = (/ -5.0, -1.0/)
    MESH_X1 = (/  5.0,  1.0/)
    BoundaryConditionsType = (/5,2,5,3/) !B wall, R transmissive, U wall, L inflow
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif

  CASE(-301) ! Shock-turbulence interaction Shu-Osher modified by Titarev and Toro
    !*Finite-volume WENO schemes for three-dimensional conservation laws, V.A. Titarev, E.F. Toro 
    TEnd    = 5.0
    Gmm     = 1.4
#ifdef MULTIFLUID
    PrimRefState4 = (/1.515695,0.523346,0.0,1.80500,0.0/)
#else
    PrimRefState4 = (/1.515695,0.523346,0.0,1.80500/)
#endif
    nElemsX = 2000 
    nElemsY = 5
    MESH_X0 = (/ -5.0, -1.0/)
    MESH_X1 = (/  5.0,  1.0/)
    BoundaryConditionsType = (/5,2,5,3/) !B wall, R transmissive, U wall, L inflow
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif


  CASE(-302) ! Shock-turbulence interaction Shu-Osher modified example 1
    TEnd    = 5.0
    Gmm     = 1.4
#ifdef MULTIFLUID
    PrimRefState4 = (/1.515695,0.523346,0.0,1.80500,0.0/)
#else
    PrimRefState4 = (/1.515695,0.523346,0.0,1.80500/)
#endif
    nElemsX = 2000 
    nElemsY = 5
    MESH_X0 = (/ -5.0, -1.0/)
    MESH_X1 = (/  5.0,  1.0/)
    BoundaryConditionsType = (/5,2,5,3/) !B wall, R transmissive, U wall, L inflow
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif

  CASE(-303) ! Shock-turbulence interaction Shu-Osher modified example 2
    TEnd    = 5.0
    Gmm     = 1.4
#ifdef MULTIFLUID
    PrimRefState4 = (/27.0/7.0,4.0/9.0*SQRT(35.0),0.0,31.0/3.0,0.0/)
#else
    PrimRefState4 = (/27.0/7.0,4.0/9.0*SQRT(35.0),0.0,31.0/3.0/)
#endif
    nElemsX = 2000 
    nElemsY = 5
    MESH_X0 = (/   -5.0, -1.0/)
    MESH_X1 = (/   15.0,  1.0/)
    BoundaryConditionsType = (/5,2,5,3/) !B wall, R transmissive, U wall, L inflow
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif

  CASE(-304) ! Shock-turbulence interaction Shu-Osher modified example 1 free BCs
    TEnd    = 5.0
    Gmm     = 1.4
    nElemsX = 2000 
    nElemsY = 5
    MESH_X0 = (/ -5.0, -1.0/)
    MESH_X1 = (/  5.0,  1.0/)
#ifdef MULTIFLUID
    PrimRefState4 = (/1.515695,0.523346,0.0,1.80500,0.0/)
#else
    PrimRefState4 = (/1.515695,0.523346,0.0,1.80500/)
#endif
    BoundaryConditionsType = (/5,2,5,2/) !B wall, R transmissive, U wall, L transmissive
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif

  CASE(-305) ! Shock-turbulence interaction Shu-Osher modified example 2 free BCs
    TEnd    = 5.0
    Gmm     = 1.4
    nElemsX = 2000 
    nElemsY = 5
    MESH_X0 = (/   -5.0, -1.0/)
    MESH_X1 = (/   15.0,  1.0/)
    BoundaryConditionsType = (/5,2,5,2/) !B wall, R transmissive, U wall, L transmissive
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif



  CASE(350) ! Woodword-Colella 
    TEnd    = 0.038
    Gmm     = 1.4
    nElemsX = 200 
    nElemsY = 5
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/1,5,1,5/) !wall
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif
  CASE(401) ! RP1 
    TEnd    = 0.2
    Gmm     = 1.4
    nElemsX = 200 
    nElemsY = 5
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/2,2,2,2/) !transmissive
    DiscontinuityLocation=0.3
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif
  CASE(402,412) ! RP2 !*412 is modified in the velocity
    TEnd    = 0.15
    Gmm     = 1.4
    nElemsX = 200 
    nElemsY = 5
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/2,2,2,2/) !transmissive
    DiscontinuityLocation=0.5
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif
  CASE(403) ! RP3
    TEnd    = 0.012
    Gmm     = 1.4
    nElemsX = 200 
    nElemsY = 5
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/2,2,2,2/) !transmissive
    DiscontinuityLocation=0.5
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif
  CASE(404) ! RP4
    TEnd    = 0.035
    Gmm     = 1.4
    nElemsX = 200 
    nElemsY = 5
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/2,2,2,2/) !transmissive
    DiscontinuityLocation=0.4
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif
  CASE(405) ! RP5
    TEnd    = 0.012
    Gmm     = 1.4
    nElemsX = 200 
    nElemsY = 5
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/2,2,2,2/) !transmissive
    DiscontinuityLocation=0.8
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif
  CASE(406) ! RP6
    TEnd    = 0.8
    Gmm     = 1.4
    nElemsX = 200 
    nElemsY = 5
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/2,2,2,2/) !transmissive
    DiscontinuityLocation=0.5
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif
  CASE(420) ! Stationary contact
    TEnd    = 2.0 !*0.5 !*0.5 !*5.0 !*50.0 !*500.0 !* 5000.0
    Gmm     = 1.4
    nElemsX = 200 
    nElemsY = 5
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/2,2,2,2/) !transmissive
    DiscontinuityLocation=0.5
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif
  CASE(421) ! Moving contact
    TEnd    = 2.0 
    Gmm     = 1.4
    nElemsX = 200 
    nElemsY = 5
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/2,2,2,2/) !transmissive
    DiscontinuityLocation=0.5
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif
  CASE(500) ! Lax
    TEnd    = 1.3
    Gmm     = 1.4
    nElemsX = 200 
    nElemsY = 5
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 10.0, 10.0/)
    Lambda_Lax = 1.0
    BoundaryConditionsType = (/2,2,2,2/) !transmissive
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif
  CASE(501) ! Lax smaller lambda
    TEnd    = 1.3
    Gmm     = 1.4
    nElemsX = 200 
    nElemsY = 5
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 10.0, 10.0/)
    Lambda_Lax = 1e-3
    BoundaryConditionsType = (/2,2,2,2/) !transmissive
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif

  !*------------------------------------------
  !*[18,28] Shock-Tube Problem - 1D in X direction Multifluid
  !*------------------------------------------
  CASE(18,28) 
    NameTest="Shock-Tube Problem - 1D in X direction Multifluid"
    TEnd    = 0.2
    Gmm     = 1.4
    nElemsX = 400
    nElemsY = 3
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/2,2,2,2/) !*Transmissive
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1=1.4
    Gmm2=1.6
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif

  !*------------------------------------------
  !*[-18] Shock-Tube Problem - 1D in X direction Multifluid reverse
  !*------------------------------------------
  CASE(-18,-28) 
    NameTest="Shock-Tube Problem - 1D in X direction Multifluid reverse"
    TEnd    = 0.2
    Gmm     = 1.4
    nElemsX = 400
    nElemsY = 3
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/2,2,2,2/) !*Transmissive
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1=1.4
    Gmm2=1.6
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif

  !*------------------------------------------
  !*[19,29] Shock-Tube Problem - 1D in Y direction Multifluid
  !*------------------------------------------
  CASE(19,29) 
    NameTest="Shock-Tube Problem - 1D in Y direction Multifluid"
    TEnd    = 0.2
    Gmm     = 1.4
    nElemsX = 120
    nElemsY = nElemsX
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/2,2,2,2/) !*Transmissive
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1=1.6
    Gmm2=1.4
    pinf1=0.0
    pinf2=0.0    
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif


  !*------------------------------------------
  !*[30] Stiff Shock-Tube Problem Multifluid
  !*------------------------------------------
  CASE(30) 
    NameTest="Stiff Shock-Tube Problem - 1D in X direction Multifluid"
    TEnd    = 0.015
    Gmm     = 1.4
    nElemsX = 400
    nElemsY = 3
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/2,2,2,2/) !*Transmissive
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1=1.4
    Gmm2=1.6
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif

  !*------------------------------------------
  !*[31] Water-Air Model Using the Stiff Equation of State Multifluid
  !*------------------------------------------
  CASE(31) 
    NameTest="Water-Air Model Using the Stiff Equation of State Multifluid"
    TEnd    = 0.00025
    Gmm     = 1.4
    nElemsX = 400
    nElemsY = 3
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/1,2,1,2/) !*Periodic BU and transmissive RL
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1=4.4
    Gmm2=1.4
    pinf1=6.0*10.0**8
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif


  !*------------------------------------------
  !*[40] Helium Bubble Multifluid
  !*------------------------------------------
  CASE(40) 
    NameTest="Helium Bubble Multifluid"
    TEnd    = 3.0
    Gmm     = 1.4
    nElemsX = 1250
    nElemsY = 500
    MESH_X0 = (/ -1.5, -0.5/)
    MESH_X1 = (/ 1.0, 0.5/)
    BoundaryConditionsType = (/5,2,5,2/) !*5 wall B U
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1=5.0/3.0
    Gmm2=1.4
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif

  !*------------------------------------------
  !*[41] R22 Bubble Multifluid
  !*------------------------------------------
  CASE(41) 
    NameTest="R22 Bubble Multifluid"
    TEnd    = 3.0
    Gmm     = 1.4
    nElemsX = 1250
    nElemsY = 500
    MESH_X0 = (/ -1.5, -0.5/)
    MESH_X1 = (/ 1.0, 0.5/)
    BoundaryConditionsType = (/5,2,5,2/) !*5 wall B U
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1=1.249
    Gmm2=1.4
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 1e-4
    NRelaxedTimesteps          = 10
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif

  !*------------------------------------------
  !*[892] Smooth periodic IC with the purpose of verifying conservation
  !*------------------------------------------
  CASE(892) 
    TEnd    = 0.5
    Gmm     = 1.4
    nElemsX = 120
    nElemsY = nElemsX
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/1.0,1.0/)
    BoundaryConditionsType = (/5,5,5,5/) !*PERIODIC BCs
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0


#ifdef RELAXATION
  !*------------------------------------------
  !*[99] Gresho vortex longer time
  !*------------------------------------------
  !*A novel full-Euler low Mach number IMEX splitting, Jonas Zeifang, Jochen Schütz, Klaus Kaiser, Andrea Beck,Maria Lukacova-Medvidova, Sebastian Noelle
  !*------------------------------------------
  CASE(99) 
    NameTest="Gresho vortex longer time"
    TEnd    = 1.0
    Gmm     = 5.0/3.0 !*Seen in many references, for example David
    nElemsX = 120
    nElemsY = nElemsX
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/1,1,1,1/) !*Periodic BCs
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0

  !*------------------------------------------
  !*[-99] Gresho vortex longer time
  !*------------------------------------------
  !*A novel full-Euler low Mach number IMEX splitting, Jonas Zeifang, Jochen Schütz, Klaus Kaiser, Andrea Beck,Maria Lukacova-Medvidova, Sebastian Noelle
  !*------------------------------------------
  CASE(-99) 
    NameTest="Gresho vortex longer time"
    TEnd    = 1.0
    Gmm     = 1.4
    nElemsX = 120
    nElemsY = nElemsX
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/1,1,1,1/) !*Periodic BCs
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0



  !*------------------------------------------
  !*[100] Gresho vortex
  !*------------------------------------------
  !*A novel full-Euler low Mach number IMEX splitting, Jonas Zeifang, Jochen Schütz, Klaus Kaiser, Andrea Beck,Maria Lukacova-Medvidova, Sebastian Noelle
  !*------------------------------------------
  CASE(100) 
    NameTest="Gresho vortex"
    TEnd    = 0.05
    Gmm     = 5.0/3.0 !*Seen in many references, for example David
    nElemsX = 120
    nElemsY = nElemsX
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/1,1,1,1/) !*Periodic BCs
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0

  !*------------------------------------------
  !*[101] Smooth Gresho vortex
  !*------------------------------------------
  !*An all speed second order IMEX relaxation scheme for the Euler equations, Andrea Thomann, Markus Zenk, Gabriella Puppo, Christian Klingenberg
  !*------------------------------------------
  CASE(101) 
    NameTest="Smooth Gresho vortex"
    TEnd    = 0.05
    Gmm     = 5.0/3.0 !*Seen in many references, for example David
    nElemsX = 120
    nElemsY = nElemsX
    MESH_X0 = (/ 0.0, 0.0/)
    MESH_X1 = (/ 1.0, 1.0/)
    BoundaryConditionsType = (/1,1,1,1/) !*Periodic BCs
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
    !*BoundaryConditionsTypeIMEXLinearSystem=(/3,3,3,3/) !*NO
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0       
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0



  !*------------------------------------------
  !*[200] 2d Low-Mach Vortex
  !*------------------------------------------
  !*A novel full-Euler low Mach number IMEX splitting, Jonas Zeifang, Jochen Schütz, Klaus Kaiser, Andrea Beck,Maria Lukacova-Medvidova, Sebastian Noelle
  !*------------------------------------------
  CASE(200) 
    NameTest="2d Low-Mach Vortex"
    TEnd    = 0.1
    Gmm     = 2.0
    nElemsX = 120
    nElemsY = nElemsX
    MESH_X0 = (/-10.0,-10.0/)
    MESH_X1 = (/10.0,10.0/)
    BoundaryConditionsType = (/1,1,1,1/) !*PERIODIC BCs
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0     
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0  

  !*------------------------------------------
  !*[201] 2d Low-Mach Vortex modified
  !*------------------------------------------
  !*A novel full-Euler low Mach number IMEX splitting, Jonas Zeifang, Jochen Schütz, Klaus Kaiser, Andrea Beck,Maria Lukacova-Medvidova, Sebastian Noelle
  !*------------------------------------------
  CASE(201) 
    NameTest="2d Low-Mach Vortex modified"
    TEnd    = 0.1
    Gmm     = 2.0
    nElemsX = 120
    nElemsY = nElemsX
    MESH_X0 = (/-10.0,-10.0/)
    MESH_X1 = (/10.0,10.0/)
    BoundaryConditionsType = (/1,1,1,1/) !*PERIODIC BCs
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0     
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0  

  !*------------------------------------------
  !*[300] Double Shear Layer
  !*------------------------------------------
  !*A novel full-Euler low Mach number IMEX splitting, Jonas Zeifang, Jochen Schütz, Klaus Kaiser, Andrea Beck,Maria Lukacova-Medvidova, Sebastian Noelle
  !*------------------------------------------
  CASE(300) 
    NameTest="Double Shear Layer"
    TEnd    = 10.0
    Gmm     = 1.4
    nElemsX = 64
    nElemsY = nElemsX
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/2.0*PI,2.0*PI/)
    BoundaryConditionsType = (/1,1,1,1/) !*PERIODIC BCs
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0     
#if defined(IMEX) || defined(IMEXMOMENTUM)
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#else
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0
#endif

  !*------------------------------------------
  !*[301] Baroclinic vorticity generation problem
  !*------------------------------------------
  !*A novel full-Euler low Mach number IMEX splitting, Jonas Zeifang, Jochen Schütz, Klaus Kaiser, Andrea Beck,Maria Lukacova-Medvidova, Sebastian Noelle
  !*------------------------------------------
  CASE(301) 
    NameTest="Baroclinic vorticity generation problem"
    TEnd    = 20.0
    Gmm     = 1.4
    nElemsX = 800
    nElemsY = 160
    EPS_LM  = 0.05
    MESH_X0 = (/-1.0/EPS_LM, 0.0/)
    MESH_X1 = (/ 1.0/EPS_LM,2.0*(1.0/EPS_LM)/5.0 /)
    BoundaryConditionsType = (/1,1,1,1/) !*PERIODIC BCs
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0     
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0  

  !*------------------------------------------
  !*[302] Kelvin-Helmholtz instability
  !*------------------------------------------
  !*All-speed numerical methods for the Euler equations via a sequential explicit time integration, Wasilij Barsukow
  !*------------------------------------------
  CASE(302) 
    NameTest="Kelvin-Helmholtz instability"
    TEnd    = 50.0
    Gmm     = 1.4
    nElemsX = 1000
    nElemsY = 500
    EPS_LM  = 0.05
    MESH_X0 = (/0.0, 0.0/)
    MESH_X1 = (/2.0, 1.0/)
    BoundaryConditionsType = (/1,1,1,1/) !*PERIODIC BCs
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
    BoundaryConditionsTypeIMEXLinearSystem=BoundaryConditionsType
#endif
#endif
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
    Gmm1 = Gmm
    Gmm2 = Gmm
    pinf1=0.0
    pinf2=0.0
#endif
#endif
    GravitationalPotentialFlag = 0     
    relaxed_dt                 = 0.0
    NRelaxedTimesteps          = 0  


#endif
#endif
  CASE DEFAULT
    ErrorMessage = "Initial condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

index=index+1
IF (nargs > index) THEN
   CALL get_command_ARGUMENT(index+1, arg)
   READ(arg, *) iarg
   nElemsX = iarg
END IF

index=index+1
IF (nargs > index) THEN
   CALL get_command_ARGUMENT(index+1, arg)
   READ(arg, *) iarg
   nElemsY = iarg
END IF


maxTimeSteps           = 100000000
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
maxLenghtJacobiCounter = 10000000
#endif
#endif

!*---------------------------------------------
!*SPACE DISCRETIZATION (RECONSTRUCTION) LEGEND
!*---------------------------------------------
!* 1=First order FV
!* 2=MUSCL
!* 3=WENO3
!* 4=WENO5
!*---------------------------------------------
!* Different MINMOD limiters
!* 20=2 = MUSCL
!* 21   = k  !*OK
!* 22   = CO !*Ok only for k=1
!* 23   = VL !*Best
!* 24   = M  !*Best
!* 25   = VA !*Not 100% coded actually
!*---------------------------------------------

Reconstruction    = -26
ReconstructionFix = Reconstruction
WhichRiemannSolver = 1                                     !* 0 Rusanov, -1 CU, 1 LDCU,   2 LDCU new antidiffusion, 3 Exact
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
WhichRiemannSolverPrimitiveSystem =1                       !* 0 Rusanov, -1 CU, 1 LDPCCU
#endif
ReconstructedVariable =2 !*0=Conserved, 1=Characteristic, 2=Primitive, -1=Characteristic Alex' and Shaoshuai


WhichSpeedEstimateForDtComputation=0 !*0=Standard, 1=Riemann


index=index+1
IF (nargs > index) THEN
  CALL get_command_ARGUMENT(index+1, arg)
  READ(arg, *) iarg
  Reconstruction = iarg
END IF

index=index+1
IF (nargs > index) THEN
  CALL get_command_ARGUMENT(index+1, arg)
  READ(arg, *) iarg
  ReconstructionFix = iarg
END IF



!*---------------------------------------------
!*TIME SCHEME LEGEND
!*---------------------------------------------
!*  1 digit
!*  1 explicit euler, 2 SSPRK2, 3 SSPRK3, 4 SSPRK64,  5 RK65
!* -1 IMEX Euler
!* 
!*  2 digits
!*  1* DeC   
!*  2* mPDeC 
!*---------------------------------------------

timescheme = 3 !*-1

index=index+1
IF (nargs > index) THEN
  CALL get_command_ARGUMENT(index+1, arg)
  READ(arg, *) iarg
  timescheme = iarg
END IF


CFL      = 0.3

index=index+1
IF (nargs > index) THEN
  CALL get_command_ARGUMENT(index+1, arg)
  READ(arg, *) iarg_real   
  CFL = iarg_real
ENDIF




WhichOutput  = 1 ! 0 Nothing, 1 Octave, 2 Tecplot, 3 Both
nOutputFiles = 1 !*1

VarNameVisu(1) = "Density"
VarNameVisu(2) = "MomentumX"
VarNameVisu(3) = "MomentumY"
VarNameVisu(4) = "Energy"
#ifdef MULTIFLUID
VarNameVisu(5) = "Level_Set_Variable" !*This is actually not going to be updated
#endif
VarNameVisu(nVar+1) = "Gravitational_Potential"

#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
VarNameVisu_W(1) = "Density"
VarNameVisu_W(2) = "VelocityX"
VarNameVisu_W(3) = "VelocityY"
VarNameVisu_W(4) = "Pressure"
VarNameVisu_W(nVar+1) = "Gravitational_Potential"
#endif

#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
VarNameVisu_W(5) = "Level_Set_Variable" !*Rho*Phi
#endif
#endif

PRINT*, "--------------------------"
PRINT*, "Test              = ", InitialCondition, TRIM(NameTest)
PRINT*, "Reconstruction    = ", Reconstruction
PRINT*, "ReconstructionFix = ", ReconstructionFix
PRINT*, "Time Scheme       = ", timescheme
SELECT CASE(ReconstructedVariable)
  CASE(0)
    PRINT*, "Reconstructed variables: Conserved"
  CASE(-1)
    PRINT*, "Reconstructed variables: Characteristic Alex' style based on interface"
  CASE(1)
    PRINT*, "Reconstructed variables: Characteristic"
  CASE(2)
    PRINT*, "Reconstructed variables: Primitive"
  CASE DEFAULT
    PRINT*, "Wrong ReconstructedVariable"
    PRINT*, ReconstructedVariable
    STOP
END SELECT
#if !defined(ACTIVEFLUX)
SELECT CASE(WhichRiemannSolver)
  CASE(0)
    PRINT*, "Riemann solver conserved system: Rusanov"
  CASE(-1)
    PRINT*, "Riemann solver conserved system: Central Upwind"
  CASE(1)
    PRINT*, "Riemann solver conserved system: Low-Dissipation Central Upwind"
  CASE(2)
    PRINT*, "Riemann solver conserved system: Low-Dissipation Central Upwind with new antidiffusion term"
  CASE(3)
    PRINT*, "Riemann solver conserved system: Exact Riemann solver"
  CASE DEFAULT
    PRINT*, "Wrong Riemann solver conserved system"
    PRINT*, WhichRiemannSolver
    STOP
END SELECT
#endif
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SELECT CASE(WhichRiemannSolverPrimitiveSystem)
  CASE(0)
    PRINT*, "Riemann solver primitive system: Rusanov"
  CASE(-1)
    PRINT*, "Riemann solver primitive system: Path Conservative Central Upwind"
  CASE(1)
    PRINT*, "Riemann solver primitive system: Low-Dissipation Path Conservative Central Upwind"
  CASE DEFAULT
    PRINT*, "Wrong Riemann solver primitive system"
    PRINT*, WhichRiemannSolverPrimitiveSystem
    STOP
END SELECT
#endif
PRINT*, "nElemsX = ", nElemsX, ", nElemsY = ", nElemsY 
PRINT*, "CFL = ", CFL
PRINT*, "--------------------------"
SELECT CASE(WhichSpeedEstimateForDtComputation)
  CASE(0)
    PRINT*, "Speed Estimate: Standard"
  CASE(1)
    PRINT*, "Speed Estimate: Riemann solver"
#ifdef RELAXATION
    PRINT*, "This speed estimate for RELAXATION is not available because the Riemann problem must be modified to keep relaxation into account"
#endif
#ifdef MULTIFLUID
    PRINT*, "This speed estimate for MULTIFLUID is not available. I need a multifluid Riemann solver."
#endif
  CASE DEFAULT
    PRINT*, "Wrong Speed estimate"
    PRINT*, WhichSpeedEstimateForDtComputation
    STOP
END SELECT
PRINT*, "--------------------------"

#ifdef SW
  PRINT*, "--------------------------"
  PRINT*, "SW setting"
  PRINT*, "p     = K*rho**gamma"
  PRINT*, "K     = ", Kappa
  PRINT*, "gamma = ", Gmm
  PRINT*, "--------------------------"
  !*Equivalence if K=g/2.0 and gamma=2.0
#endif

#ifdef PATANKAR
  PRINT*, "--------------------------"
  PRINT*, "Patankar active"
  PRINT*, "--------------------------"
#endif

#ifdef WELLBALANCED
  PRINT*, "--------------------------"
  PRINT*, "Well-balanced active"
  PRINT*, "--------------------------"
#endif


#ifdef CENTEREDPRIMITIVE
  PRINT*, "--------------------------"
  PRINT*, "Centered primitive active"
#ifdef RECONSTRUCTFROMCONSERVED
  PRINT*, "--------------------------"
  PRINT*, "Reconstruct interface values from conserved variables"
#else
  PRINT*, "--------------------------"
  PRINT*, "Reconstruct interface values from primitive variables"
#endif
#ifdef MOMENTUMINPRIMITIVEVARIABLES
  PRINT*, "--------------------------"
  PRINT*, "Momentum in primitive variables"
#endif
#ifdef FIRSTORDERPRIMITIVE
  PRINT*, "--------------------------"
  PRINT*, "First order primitive (cheaper shock detector)"
#endif
#ifdef EVOLVERHOWITHSAMESCHEME
  PRINT*, "--------------------------"
  PRINT*, "Evolve rho with same scheme"
#endif
#ifdef PATHCONSERVATIVESHOCKDETECTION
  K_coefficient_PCSD = 0.025
  index=index+1
  IF (nargs > index) THEN
    CALL get_command_ARGUMENT(index+1, arg)
    READ(arg, *) iarg_real
    K_coefficient_PCSD = iarg_real
  END IF
  PRINT*, "--------------------------"
  PRINT*, "Shock detection based on Path Conservative"
  PRINT*, "...with K=", K_coefficient_PCSD
  WhichVariableForLimiting_PCSD = 3
  index=index+1
  IF (nargs > index) THEN
    CALL get_command_ARGUMENT(index+1, arg)
    READ(arg, *) iarg
    WhichVariableForLimiting_PCSD = iarg
  END IF
  SELECT CASE(WhichVariableForLimiting_PCSD)
    CASE(1)
      PRINT*, "...limiting on density"
    CASE(2)
      PRINT*, "...limiting on momentum"
    CASE(3)
      PRINT*, "...limiting on energy"
    CASE DEFAULT
      PRINT*, "Wrong choice for limiting variable"
      PRINT*, WhichVariableForLimiting_PCSD
      STOP
  END SELECT
#ifdef SMOOTHINGLIMITER
      PRINT*, "...smoothing limiter active"
#else
      PRINT*, "...smoothing limiter NOT active"
#endif
#ifdef DETECTEACHSTAGE
      PRINT*, "...detect at each stage"
#else
      PRINT*, "...detect only at the beginning of the time step"
#endif
#endif
#ifdef DONOTOVERWRITEPRIMITIVE
      PRINT*, "Not overwriting primitive, basically for debugging"
#else
      PRINT*, "Overwriting primitive as it should be"
#endif
  PRINT*, "--------------------------"
#endif




#ifdef ACTIVEFLUX
  PRINT*, "--------------------------"
  PRINT*, "Active flux active"

  !*---------------------------------------------
  !*POSTPROCESSING LEGEND
  !*---------------------------------------------
  !*  0 nothing
  !*  1 minmod, first x and then y
  !*  2 minmod, first y and then x
  !*  The following depend on eps
  !*  3 minmod, first x and then y WORKING BUT NOT CONSERVATIVE
  !*  4 minmod, first y and then x WORKING BUT NOT CONSERVATIVE
  !*  5 minmod, first x and then y CONSERVATIVE
  !*  6 minmod, first y and then x CONSERVATIVE
  !*---------------------------------------------

  AF_PostProcessing = 1
  index=index+1
  IF (nargs > index) THEN
    CALL get_command_ARGUMENT(index+1, arg)
    READ(arg, *) iarg
    AF_PostProcessing = iarg
  END IF

#ifdef PRIMITIVEONLY
  PRINT*, "Evolution of primitive variables only"
#else
  PRINT*, "Postprocessing = ", AF_PostProcessing
#if (defined(OTHERPOSTPROCESSING) + defined(FOURINPUTSPOSTPROCESSING) + defined(ALINAPOSTPROCESSING) + defined(SIDEDPOSTPROCESSING)) > 1
  PRINT*, "Impossible to have multiple post-processing flags active"
  STOP
#endif


#endif

  PRINT*, "--------------------------"
#endif

#ifdef PRIMITIVEONLY
#ifndef ACTIVEFLUX
  PRINT*, "Impossible to evolve only primitive variables without Active Flux flag."
  STOP
#endif
#endif

#ifdef RELAXATION
  PRINT*, "--------------------------"
  PRINT*, "Relaxation"


  EPS_LM=1.0
  index=index+1
  IF (nargs > index) THEN
    CALL get_command_ARGUMENT(index+1, arg)
    READ(arg, *) iarg_real   
    EPS_LM = iarg_real
  ENDIF


  PRINT*, "eps = ", EPS_LM

  IF ((InitialCondition .EQ. 301) .AND. (EPS_LM .NE. 0.05) ) THEN
    PRINT*, "For this test eps mus be 0.05"
    STOP
  ENDIF

  PRINT*, "--------------------------"
#endif


#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#ifdef IMEX
  SAFETY_K_TO_AVOID_HYPERBOLICITY_LOSS=0.0
  index=index+1
  IF (nargs > index) THEN
    CALL get_command_ARGUMENT(index+1, arg)
    READ(arg, *) iarg_real   
    SAFETY_K_TO_AVOID_HYPERBOLICITY_LOSS = iarg_real
  ENDIF
#endif
#endif

#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#ifdef IMEX
  PRINT*, "--------------------------"
  PRINT*, "IMEX active"
#ifdef NOMODIFICATIONWAVESPEEDIMEX
  PRINT*, "...without hyperbolicity trick. K=0, i.e., no modification of the wave speeds in IMEX"
  PRINT*, "YOU'RE NOT SUPPOSED TO USE THIS"
  STOP
#elif HYPERBOLICITYTRICKMODIFIED
  PRINT*, "...with hyperbolicity trick modified, K is a function of Mach, i.e., epsilon"
#else
  PRINT*, "...with original hyperbolicity trick"
  PRINT*, "...with safety K to avoid hyperbolicity loss", SAFETY_K_TO_AVOID_HYPERBOLICITY_LOSS
  PRINT*, "YOU'RE NOT SUPPOSED TO USE THIS"
  STOP
#endif
#ifdef IMEXL2FULLYUPWINDED
  PRINT*, "IMEX DeC with L^2 fully upwinded..."
  IF (timescheme .EQ. -1) THEN
    PRINT*, "...This feature only works with DeC"
    STOP
  END IF
  PRINT*, "ACTUALLY, THIS FEATURE NEVER WORKS PROPERLY"
  STOP
#endif
#ifdef ALTERNATIVEFORMULATIONPRESSURE
  PRINT*, "ALTERNATIVEFORMULATIONPRESSURE not compatible with IMEX splitting"
  STOP
#endif
#ifdef BLENDINGDT
  PRINT*, "...with blending of dt implicit and implicit"
  PRINT*, "YOU'RE NOT SUPPOSED TO USE THIS"
  STOP
#endif
  PRINT*, "--------------------------"
#endif

#ifdef IMEXMOMENTUM
  PRINT*, "--------------------------"
  PRINT*, "IMEXMOMENTUM active"
#ifndef MOMENTUMINPRIMITIVEVARIABLES
  PRINT*, "IMEXMOMENTUM requires MOMENTUMINPRIMITIVEVARIABLES"
  STOP
#endif
#ifndef PRESSUREFORMULATIONFORIMEXMOMENTUM
  PRINT*, "IMEXMOMENTUM requires PRESSUREFORMULATIONFORIMEXMOMENTUM"
  STOP
#endif
  PRINT*, "--------------------------"
#endif
#endif


#ifdef ACTIVEFLUX
#ifdef MOMENTUMINPRIMITIVEVARIABLES
  PRINT*, "--------------------------"
  PRINT*, "Momentum in primitive variables"
  PRINT*, "--------------------------"
#endif
#ifdef RECONSTRUCTIONVELOCITY
  PRINT*, "--------------------------"
  PRINT*, "Reconstruction velocity in primitive variables"
  PRINT*, "--------------------------"
#ifndef MOMENTUMINPRIMITIVEVARIABLES
  PRINT*, "RECONSTRUCTIONVELOCITY is kind of redundant and useless if I do not have MOMENTUMINPRIMITIVEVARIABLES"
  PRINT*, "Disable it"
  STOP
#endif
#endif
#ifdef ALTERNATIVEFORMULATIONPRESSURE
  PRINT*, "--------------------------"
  PRINT*, "Alternative formulation for pressure equation, the one which is not good for IMEX splitting"
  PRINT*, "--------------------------"
#endif
#ifdef PRESSUREFORMULATIONFORIMEXMOMENTUM
  PRINT*, "--------------------------"
  PRINT*, "Alternative formulation for pressure equation, the one which is suitable for IMEX splitting and momentum equation"
  PRINT*, "--------------------------"
#endif
#if (defined(ALTERNATIVEFORMULATIONPRESSURE) + defined(PRESSUREFORMULATIONFORIMEXMOMENTUM) ) > 1
  PRINT*, "Impossible to have multiple PRESSUREFORMULATION flags active"
  STOP
#endif
#endif


SELECT CASE (timescheme)
  CASE(11,12,21,22,-2,-12,-22,-32,-52,-82,-102)   
    MstepsMax=2 !*Order 1,2
  CASE(13,14,23,24,-3,-4,-13,-14,-23,-33,-24,-34,-103,-104) 
    MstepsMax=3 !*Order 3,4
  CASE(15,25,-5,-15,-25,-35,-105) 
    MstepsMax=4 !*Order 5
  CASE DEFAULT
    MstepsMax=1
END SELECT 


#ifdef ACTIVEFLUX

#if defined(PRIMITIVEONLY) && defined(IMPLICITUPDATECONSERVED)
PRINT*, "PRIMITIVEONLY incompatible with IMPLICITUPDATECONSERVED"
STOP
#endif

#ifdef IMEX

#ifdef NORMALIZEALL
PRINT*, "Normalization of linear systems by all"
#elif defined(NORMALIZELINEARSYSTEMS)
PRINT*, "Normalization of linear systems by dx only"
PRINT*, "I prefer to normalize by all, so I do not let it run in this mode"
STOP
#else
PRINT*, "No normalization"
PRINT*, "I prefer to normalize by all, so I do not let it run in this mode"
STOP
#endif
PRINT*, "--------------------------"


#if defined(NORMALIZEALL) && defined(NORMALIZELINEARSYSTEMS)
  PRINT*, "Double normalization not possible"
  PRINT*, "Go to makefile and get rid of one of the two"
  STOP
#endif

#endif

#endif


#ifdef MULTIFLUID
PRINT*, "--------------------------"
PRINT*, "Multifluid active"
PRINT*, "Gamma 1", Gmm1
PRINT*, "Gamma 2", Gmm2
#ifndef ACTIVEFLUX
PRINT*, "Impossible to run Multifluid without ACTIVEFLUX"
STOP
#endif
PRINT*, "...with N troubled neighbours: ", N_Troubled_Neighbors
#if defined(RELAXATION) || defined(IMEX)
PRINT*, "Multifluid incompatible with RELAXATION or IMEX"
STOP
#endif

#ifdef LEVELSETINU
PRINT*, "...evolving also the level set variable in U"
#endif

#ifdef PWCINTROUBLEDCELLS
PRINT*, "...piecewise constant reconstruction in troubled cells"
#endif

#ifdef RECONSTRUCTIONPHI
PRINT*, "...reconstruction of phi rather than rophi"
IF ( ReconstructedVariable .EQ. -1 ) THEN
    PRINT*, "RECONSTRUCTIONPHI flag incompatible with reconstruction of characteristic variables from Alex and Shaoshuai"
    STOP
END IF
#endif

#if defined(MOMENTUMINPRIMITIVEVARIABLES) && !defined(RECONSTRUCTIONVELOCITY)
IF ( ReconstructedVariable .EQ. -1 ) THEN
    PRINT*, "Reconstruction of Characteristic variables not set up for MOMENTUMINPRIMITIVEVARIABLES without RECONSTRUCTIONVELOCITY"
    STOP
END IF
#endif



#ifdef PRESSUREFORMULATIONFORIMEXMOMENTUM
PRINT*, "PRESSUREFORMULATIONFORIMEXMOMENTUM not yet coded for MULTIFLUID"
STOP
#endif

PRINT*, "--------------------------"

#endif

#if defined(CENTEREDPRIMITIVE) &&  defined(ACTIVEFLUX)
PRINT*, "CENTEREDPRIMITIVE and ACTIVEFLUX are incompatible"
STOP
#endif

#if defined(CENTEREDPRIMITIVE) &&  defined(RECONSTRUCTIONVELOCITY)
PRINT*, "CENTEREDPRIMITIVE and RECONSTRUCTIONVELOCITY are incompatible"
STOP
#endif

#if defined(FIRSTORDERPRIMITIVE) && !defined(RECONSTRUCTFROMCONSERVED)
PRINT*, "FIRSTORDERPRIMITIVE requires RECONSTRUCTFROMCONSERVED"
STOP
#endif

#if defined(EVOLVERHOWITHSAMESCHEME) && !defined(CENTEREDPRIMITIVE)
PRINT*, "EVOLVERHOWITHSAMESCHEME requires CENTEREDPRIMITIVE"
STOP
#endif


#if defined(EVOLVERHOWITHSAMESCHEME) && defined(CENTEREDPRIMITIVE)
  IF (((WhichRiemannSolverPrimitiveSystem .NE. -1) .AND. (WhichRiemannSolver .NE. -1) ) &
    & .AND. ( (WhichRiemannSolverPrimitiveSystem .NE. 1) .AND. ((WhichRiemannSolver .NE. 1) .OR. (WhichRiemannSolver .NE. 2)) ) ) THEN
      PRINT*, "EVOLVERHOWITHSAMESCHEME requires either"
      PRINT*, "CU for both primitive and conserved"
      PRINT*, "or"
      PRINT*, "LDCU for both primitive and conserved"
  END IF
#endif




#if defined(EVOLVERHOWITHSAMESCHEME) && defined(OLDEVOLVERHOWITHSAMESCHEME)
PRINT*, "EVOLVERHOWITHSAMESCHEME and OLDEVOLVERHOWITHSAMESCHEME are incompatible"
PRINT*, "Select one of the two indicators."
STOP 
#endif

#ifdef MULTIFLUID
IF (ReconstructedVariable .EQ. 1) THEN
#ifndef RECONSTRUCTIONPHI
  PRINT*, "Reconstruction of characteristic variables with MULTIFLUID requires RECONSTRUCTIONPHI"
  STOP
#endif
#ifdef MOMENTUMINPRIMITIVEVARIABLES
#ifndef RECONSTRUCTIONVELOCITY
  PRINT*, "Reconstruction of characteristic variables with MULTIFLUID and MOMENTUMINPRIMITIVEVARIABLES requires RECONSTRUCTIONVELOCITY"
  STOP
#endif
#endif
ELSE IF (ReconstructedVariable .EQ. 0) THEN
  PRINT*, "You may have selected by chance reconstruction of conserved variables with multifluid"
  PRINT*, "I'm stopping before running."
  STOP
END IF
#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE InitializeParameters
!===============================================================================!
!
!
!
!===============================================================================!
END MODULE MOD_Parameters
!-------------------------------------------------------------------------------!
