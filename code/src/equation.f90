!===============================================================================!
MODULE MOD_Equation
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE ExactFunction
  MODULE PROCEDURE ExactFunction
END INTERFACE

INTERFACE ExactFunctionWB
  MODULE PROCEDURE ExactFunctionWB
END INTERFACE

INTERFACE SourceTerms
  MODULE PROCEDURE SourceTerms
END INTERFACE

INTERFACE SourceTerm_Conserved_INPUT_CONSERVED
  MODULE PROCEDURE SourceTerm_Conserved_INPUT_CONSERVED
END INTERFACE

#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
INTERFACE SourceTerm_Primitive_INPUT_PRIMITIVE
  MODULE PROCEDURE SourceTerm_Primitive_INPUT_PRIMITIVE
END INTERFACE
#endif

INTERFACE BoundaryConditions
  MODULE PROCEDURE BoundaryConditions
END INTERFACE

INTERFACE TimeStep
  MODULE PROCEDURE TimeStep
END INTERFACE

INTERFACE RiemannSolver
  MODULE PROCEDURE RiemannSolver
END INTERFACE

INTERFACE EvaluateFlux1D
  MODULE PROCEDURE EvaluateFlux1D
END INTERFACE

INTERFACE ConsToPrim
  MODULE PROCEDURE ConsToPrim
END INTERFACE

INTERFACE PrimToCons
  MODULE PROCEDURE PrimToCons
END INTERFACE

INTERFACE Gravitational_Potential   
  MODULE PROCEDURE Gravitational_Potential   
END INTERFACE

#ifndef MULTIFLUID
INTERFACE LeftEigenvectors
  MODULE PROCEDURE LeftEigenvectors
END INTERFACE

INTERFACE RightEigenvectors
  MODULE PROCEDURE RightEigenvectors
END INTERFACE

INTERFACE LeftEigenvectorsObtainedFromRotationalInvariance
  MODULE PROCEDURE LeftEigenvectorsObtainedFromRotationalInvariance
END INTERFACE

INTERFACE RightEigenvectorsObtainedFromRotationalInvariance
  MODULE PROCEDURE RightEigenvectorsObtainedFromRotationalInvariance
END INTERFACE


INTERFACE LeftEigenvectorsXDirection
  MODULE PROCEDURE LeftEigenvectorsXDirection
END INTERFACE

INTERFACE RightEigenvectorsXDirection
  MODULE PROCEDURE RightEigenvectorsXDirection
END INTERFACE

INTERFACE LeftEigenvectorsYDirection
  MODULE PROCEDURE LeftEigenvectorsYDirection
END INTERFACE

INTERFACE RightEigenvectorsYDirection
  MODULE PROCEDURE RightEigenvectorsYDirection
END INTERFACE
#endif

INTERFACE LeftEigenvectorsRotationalInvariancePrimitiveSystem
  MODULE PROCEDURE LeftEigenvectorsRotationalInvariancePrimitiveSystem
END INTERFACE

INTERFACE RightEigenvectorsRotationalInvariancePrimitiveSystem
  MODULE PROCEDURE RightEigenvectorsRotationalInvariancePrimitiveSystem
END INTERFACE


INTERFACE LeftEigenvectorsXDirectionPrimitiveSystem
  MODULE PROCEDURE LeftEigenvectorsXDirectionPrimitiveSystem
END INTERFACE

INTERFACE RightEigenvectorsXDirectionPrimitiveSystem
  MODULE PROCEDURE RightEigenvectorsXDirectionPrimitiveSystem
END INTERFACE

INTERFACE LeftEigenvectorsYDirectionPrimitiveSystem
  MODULE PROCEDURE LeftEigenvectorsYDirectionPrimitiveSystem
END INTERFACE

INTERFACE RightEigenvectorsYDirectionPrimitiveSystem
  MODULE PROCEDURE RightEigenvectorsYDirectionPrimitiveSystem
END INTERFACE


INTERFACE TransitionMatrixConsToPrim
  MODULE PROCEDURE TransitionMatrixConsToPrim
END INTERFACE

INTERFACE TransitionMatrixPrimToCons
  MODULE PROCEDURE TransitionMatrixPrimToCons
END INTERFACE

#if defined(CENTEREDPRIMITIVE) && defined(PATHCONSERVATIVESHOCKDETECTION)
INTERFACE Impose_BC_on_Troubled_Cell_INPUT
  MODULE PROCEDURE Impose_BC_on_Troubled_Cell_INPUT
END INTERFACE
#endif

#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
INTERFACE NormalFlux1D
  MODULE PROCEDURE NormalFlux1D
END INTERFACE

INTERFACE NumericalFluxPrimitiveSystem
  MODULE PROCEDURE NumericalFluxPrimitiveSystem
END INTERFACE

INTERFACE PathConservativeSurfaceContribution
  MODULE PROCEDURE PathConservativeSurfaceContribution
END INTERFACE

INTERFACE BoundaryConditions_On_Cell_Centered_Data_Primitive_INPUT
  MODULE PROCEDURE BoundaryConditions_On_Cell_Centered_Data_Primitive_INPUT
END INTERFACE


#ifdef IMEX
INTERFACE Compute_min_p_max_ro
  MODULE PROCEDURE Compute_min_p_max_ro
END INTERFACE
#endif
#endif

#ifdef CENTEREDPRIMITIVE
INTERFACE Impose_BC_on_WCt
  MODULE PROCEDURE Impose_BC_on_WCt
END INTERFACE

INTERFACE Impose_BC_on_WC
  MODULE PROCEDURE Impose_BC_on_WC
END INTERFACE
#endif

#ifdef ACTIVEFLUX
INTERFACE Impose_BC_on_Wt
  MODULE PROCEDURE Impose_BC_on_Wt
END INTERFACE

INTERFACE Impose_BC_on_W
  MODULE PROCEDURE Impose_BC_on_W
END INTERFACE

#ifdef MULTIFLUID
INTERFACE Get_Gamma   
  MODULE PROCEDURE Get_Gamma   
END INTERFACE

INTERFACE Get_Pressure_Infinity   
  MODULE PROCEDURE Get_Pressure_Infinity
END INTERFACE

INTERFACE Impose_BC_on_Troubled_Cell_U   
  MODULE PROCEDURE Impose_BC_on_Troubled_Cell_U   
END INTERFACE

INTERFACE Impose_BC_on_Troubled_Cell_W   
  MODULE PROCEDURE Impose_BC_on_Troubled_Cell_W   
END INTERFACE

INTERFACE Impose_BC_on_Fluid_Cell_U   
  MODULE PROCEDURE Impose_BC_on_Fluid_Cell_U   
END INTERFACE
#endif

#endif

INTERFACE Switching_Function
  MODULE PROCEDURE Switching_Function
END INTERFACE

!-------------------------------------------------------------------------------!
PUBLIC :: ExactFunction
PUBLIC :: ExactFunctionWB
PUBLIC :: SourceTerms
PUBLIC :: SourceTerm_Conserved_INPUT_CONSERVED
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
PUBLIC :: SourceTerm_Primitive_INPUT_PRIMITIVE
#endif
PUBLIC :: BoundaryConditions
PUBLIC :: TimeStep
PUBLIC :: RiemannSolver
PUBLIC :: EvaluateFlux1D
PUBLIC :: ConsToPrim
PUBLIC :: PrimToCons
PUBLIC :: Gravitational_Potential
#ifndef MULTIFLUID
PUBLIC :: LeftEigenvectors
PUBLIC :: RightEigenvectors
PUBLIC :: LeftEigenvectorsObtainedFromRotationalInvariance
PUBLIC :: RightEigenvectorsObtainedFromRotationalInvariance
PUBLIC :: LeftEigenvectorsXDirection
PUBLIC :: RightEigenvectorsXDirection
PUBLIC :: LeftEigenvectorsYDirection
PUBLIC :: RightEigenvectorsYDirection
#endif
PUBLIC :: LeftEigenvectorsRotationalInvariancePrimitiveSystem
PUBLIC :: RightEigenvectorsRotationalInvariancePrimitiveSystem
PUBLIC :: LeftEigenvectorsXDirectionPrimitiveSystem
PUBLIC :: RightEigenvectorsXDirectionPrimitiveSystem
PUBLIC :: LeftEigenvectorsYDirectionPrimitiveSystem
PUBLIC :: RightEigenvectorsYDirectionPrimitiveSystem
PUBLIC :: TransitionMatrixConsToPrim
PUBLIC :: TransitionMatrixPrimToCons
PUBLIC :: BoundaryConditions_On_Cell_Centered_Data_Primitive_INPUT

#if defined(CENTEREDPRIMITIVE) && defined(PATHCONSERVATIVESHOCKDETECTION)
PUBLIC :: Impose_BC_on_Troubled_Cell_INPUT
#endif

#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
PUBLIC :: NormalFlux1D
PUBLIC :: NumericalFluxPrimitiveSystem
PUBLIC :: PathConservativeSurfaceContribution

#ifdef IMEX
PUBLIC :: Compute_min_p_max_ro
#endif
#endif

#ifdef ACTIVEFLUX
PUBLIC :: Impose_BC_on_Wt
PUBLIC :: Impose_BC_on_W

#ifdef MULTIFLUID
PUBLIC :: Get_Gamma
PUBLIC :: Get_Pressure_Infinity
PUBLIC :: Impose_BC_on_Troubled_Cell_U
PUBLIC :: Impose_BC_on_Troubled_Cell_W
PUBLIC :: Impose_BC_on_Fluid_Cell_U
#endif
#endif

#ifdef CENTEREDPRIMITIVE
PUBLIC :: Impose_BC_on_WCt
PUBLIC :: Impose_BC_on_WC
#endif

PUBLIC :: Switching_Function
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
SUBROUTINE ExactFunction(WhichInitialCondition,t,x,Cons)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: PI
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: MESH_SX
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DENSITY
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#ifdef SW
USE MOD_FiniteVolume2D_vars,ONLY: Gravity
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
#endif
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
USE MOD_FiniteVolume2D_vars,ONLY: DiscontinuityLocation
USE MOD_FiniteVolume2D_vars,ONLY: Lambda_Lax
USE MOD_FiniteVolume2D_vars,ONLY: Ma_shock
USE MOD_FiniteVolume2D_vars,ONLY: Ma_vortex
USE MOD_FiniteVolume2D_vars,ONLY: R_EOS
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
INTEGER,INTENT(IN) :: WhichInitialCondition
REAL,INTENT(IN)    :: t
REAL,INTENT(IN)    :: x(1:nDims)
REAL,INTENT(OUT)   :: Cons(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: Prim(1:nVar)
REAL               :: xc(2), xm(2), r, r0, hl, hr, r2, r20
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!
!*OUR VARIABLES
REAL               :: Omega, Jamma, u_inf, v_inf, h_inf, DeltaH
INTEGER            :: power
REAL               :: ro_inf, p_inf, beta, delta_u, delta_v, delta_T, delta_ro
REAL               :: xmxc(1:2), x_wrt_BL(1:2), x_wrt_BL_bm(1:2), x_0(1:2), x_d(1:2)
REAL               :: ri
REAL               :: p0, sintheta, costheta, ro_bar, vmax, energy, delta_energy, theta
REAL               :: u0, ro0
REAL               :: level_set
REAL               :: a, b, vm, vtheta
REAL               :: ro1, u1, v1, p1, T1
REAL               :: ro2, u2, v2, p2
REAL               :: Acoeff, Bcoeff, Ccoeff
REAL               :: Temp
REAL               :: ss, ro, p

Cons = 0.0
Prim = 0.0
SELECT CASE(WhichInitialCondition)
#ifdef SW
  !*------------------------------------------
  !*[1] Unsteady smooth vortex for SW
  !*------------------------------------------
  CASE(1) 
    u_inf = 2.
    v_inf = 3.
    H_inf=1.
    r0 = 1.

    xm(1) = MESH_X0(1)+0.5*MESH_SX(1)
    xm(2) = MESH_X0(2)+0.5*MESH_SX(2)
    xc(1) = MODULO( x(1)-u_inf*t-MESH_X0(1) , Mesh_SX(1) ) + MESH_X0(1)-xm(1)
    xc(2) = MODULO( x(2)-v_inf*t-MESH_X0(2) , Mesh_SX(2) ) + MESH_X0(2)-xm(2)
    r     = (xc(1)**2 + xc(2)**2)
    Omega = sqrt(2.*Gravity*hDerivSmoothAuxiliary(r))
    
    Prim(1) = H_inf
    Prim(2) = u_inf
    Prim(3) = v_inf

    IF (r .LT. 1) THEN
      Prim(1) = hSmoothAuxiliary(r)
      Prim(2) = Prim(2)+Omega*(+xc(2))
      Prim(3) = Prim(3)+Omega*(-xc(1))
    END IF

    Prim(4)= Kappa*Prim(1)**Gmm

#ifdef MULTIFLUID
    Prim(5)=0.0
#endif


#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)
#else
  !*------------------------------------------
  !*[2] Steady isentropic vortex
  !*------------------------------------------
  CASE(2) 

    u_inf=0.0
    v_inf=0.0

    !*Center of the vortex
    xc=0.5*(MESH_X1+MESH_X0)

    !*Coordinates from the center of the vortex
    xmxc=x-xc

    !*Distance squared from the center of the vortex
    r2=xmxc(1)**2+xmxc(2)**2

    !*Vortex amplitude
    beta=5.0 !*5.0 0.1

    delta_u=beta/(2.0*PI)*EXP( 0.5*( 1.0-r2 ) )*( -xmxc(2) )
    delta_v=beta/(2.0*PI)*EXP( 0.5*( 1.0-r2 ) )*xmxc(1)
    delta_T=-(Gmm-1.0)*beta**2/(8.0*Gmm*Pi**2)*EXP( 1.0-r2 )

    Prim(1)=(1.0+delta_T)**( 1.0 / (Gmm-1.0) )
    Prim(2)=delta_u
    Prim(3)=delta_v
    Prim(4)=(1.0+delta_T)**( Gmm / (Gmm-1.0) )

#ifdef MULTIFLUID
    Prim(5)=0.0
#endif

#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)

  !*------------------------------------------
  !*[3] Unsteady isentropic vortex
  !*------------------------------------------
  CASE(3) 

    u_inf=1.0
    v_inf=1.0

    !*Original center of the vortex before the moevement
    xc=0.5*(MESH_X1+MESH_X0)

    !*We want to get the initial position of x before the movement

    !*x with respect to bottom-left corner
    x_wrt_BL=x-MESH_X0

    !*x with respect to bottom-left corner before movement
    x_wrt_BL_bm(1)=x_wrt_BL(1)-u_inf*t
    x_wrt_BL_bm(2)=x_wrt_BL(2)-v_inf*t

    !*This is the position before movement modulo the length of the domain
    !*NB: MODULO RESULT IS ALWAYS POSITIVE
    x_wrt_BL_bm(1)=MODULO( x_wrt_BL_bm(1), MESH_SX(1) )
    x_wrt_BL_bm(2)=MODULO( x_wrt_BL_bm(2), MESH_SX(2) )

    !*This is the initial position
    x_0=MESH_X0+x_wrt_BL_bm

    !*Distance squared from the center of the vortex at the initial time
    x_d=x_0-xc

    r2=x_d(1)**2+x_d(2)**2

    !*Vortex amplitude
    beta=5.0 !*5.0 0.1

    delta_u=beta/(2.0*PI)*EXP( 0.5*( 1.0-r2 ) )*( -x_d(2) )
    delta_v=beta/(2.0*PI)*EXP( 0.5*( 1.0-r2 ) )*x_d(1)
    delta_T=-(Gmm-1.0)*beta**2/(8.0*Gmm*Pi**2)*EXP( 1.0-r2 )

    Prim(1)=(1.0+delta_T)**( 1.0 / (Gmm-1.0) )
    Prim(2)=u_inf+delta_u
    Prim(3)=v_inf+delta_v
    Prim(4)=(1.0+delta_T)**( Gmm / (Gmm-1.0) )

#ifdef MULTIFLUID
    Prim(5)=0.0
#endif

#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)

  !*------------------------------------------
  !*[4] Advection of smooth density
  !*------------------------------------------
  CASE(4)
    u_inf = 1.0
    v_inf = -0.5
    p_inf = 1.0
    Prim(1)=1.0+0.5*SIN( 4.0*Pi*( x(1)+x(2)-t*(u_inf+v_inf) ) )
    Prim(2)=u_inf
    Prim(3)=v_inf
    Prim(4)=p_inf

#ifdef MULTIFLUID
    Prim(5)=0.0
#endif

#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)

  !*------------------------------------------
  !*[5] 1D advection of smooth density in X direction
  !*------------------------------------------
  CASE(5)
    u_inf = 1.0
    v_inf = 0.0 !*<=== NB: 0 
    p_inf = 1.0
    Prim(1)=1.0+0.5*SIN( 4.0*Pi*( x(1)-t*u_inf ) ) !*<=== NB: Independent on y
    Prim(2)=u_inf
    Prim(3)=v_inf
    Prim(4)=p_inf

#ifdef MULTIFLUID
    Prim(5)=0.0
#endif

#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)


  CASE(151) ! sin4 wave
    !*=================================
    !*TVD Fluxes for the High-Order ADER Schemes, Toro, Titarev
    !*=================================
    u_inf = 1.0
    v_inf = 0.0 !*<=== NB: 0 
    p_inf = 1.0
    Prim(1)=2.0+SIN(PI*(x(1)-v_inf*t))**4 !*<=== NB: Independent on y
    Prim(2)=u_inf
    Prim(3)=v_inf
    Prim(4)=p_inf

#ifdef MULTIFLUID
    Prim(5)=0.0
#endif

#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)



  !*------------------------------------------
  !*[6] 1D advection of smooth density in Y direction
  !*------------------------------------------
  CASE(6)
    u_inf = 0.0 !*<=== NB: 0
    v_inf = 1.0 
    p_inf = 1.0
    Prim(1)=1.0+0.5*SIN( 4.0*Pi*( x(2)-t*v_inf ) ) !*<=== NB: Independent on x
    Prim(2)=u_inf
    Prim(3)=v_inf
    Prim(4)=p_inf

#ifdef MULTIFLUID
    Prim(5)=0.0
#endif

#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)


  !*------------------------------------------
  !*[7] Sod - Explosion problem
  !*------------------------------------------
  CASE(7,-7,-70,77) 

    r  = SQRT(x(1)**2 + x(2)**2) !*Radial distance from the center
    ri = 0.4

    Prim(2)=0.
    Prim(3)=0.
    Prim(4)=1.

    IF (r .LT. ri) THEN
      Prim(1)=1.0
      Prim(4)=1.0
    ELSE
      Prim(1)=0.125
      Prim(4)=0.1

    END IF

#ifdef MULTIFLUID
    Prim(5)=0.0
#endif

#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)

  !*------------------------------------------
  !*[8] Sod - 1D in X direction
  !*------------------------------------------
  CASE(8,-8,-80,88) 


    Prim(2)=0.
    Prim(3)=0.
    Prim(4)=1.

    IF (x(1) .LT. 0.5) THEN
      Prim(1)=1.0
      Prim(4)=1.0
    ELSE
      Prim(1)=0.125
      Prim(4)=0.1

    END IF

#ifdef MULTIFLUID
    Prim(5)=0.0
#endif

#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)

  !*------------------------------------------
  !*[9] Sod - 1D in Y direction
  !*------------------------------------------
  CASE(9) 


    Prim(2)=0.
    Prim(3)=0.
    Prim(4)=1.

    IF (x(2) .LT. 0.5) THEN
      Prim(1)=1.0
      Prim(4)=1.0
    ELSE
      Prim(1)=0.125
      Prim(4)=0.1

    END IF

#ifdef MULTIFLUID
    Prim(5)=0.0
#endif

#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)

  !*------------------------------------------
  !*[10] Implosion problem
  !*------------------------------------------
  CASE(10) 

    r=ABS(x(1))+ABS(x(2))

    Prim(2)=0.
    Prim(3)=0.



    IF (r .LT. 0.15) THEN
      Prim(1)=0.125
      Prim(4)=0.14
    ELSE
      Prim(1)=1.0
      Prim(4)=1.0
    END IF

#ifdef MULTIFLUID
    Prim(5)=0.0
#endif

#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)


  !*------------------------------------------
  !*[11] 2D-Riemann problem
  !*------------------------------------------
  CASE(11) 

    IF ((x(1) .GE. 1.0) .AND. (x(2) .GE. 1.0)) THEN
      Prim(1)=1.5
      Prim(2)=0.0
      Prim(3)=0.0
      Prim(4)=1.5
    ELSE IF ((x(1) .LE. 1.0) .AND. (x(2) .GE. 1.0)) THEN
      Prim(1)=0.5323
      Prim(2)=1.206
      Prim(3)=0.0
      Prim(4)=0.3
    ELSE IF ((x(1) .LE. 1.0) .AND. (x(2) .LE. 1.0)) THEN
      Prim(1)=0.138
      Prim(2)=1.206
      Prim(3)=1.206
      Prim(4)=0.029
    ELSE
      Prim(1)=0.5323
      Prim(2)=0.0
      Prim(3)=1.206
      Prim(4)=0.3
    END IF

#ifdef MULTIFLUID
    Prim(5)=0.0
#endif

#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)

  !*------------------------------------------
  !*[12] Rayleigh-Taylor instability
  !*------------------------------------------
  CASE(12) 

    IF (x(2) .LT. 0.5) THEN
      Prim(1)=2.0
      Prim(4)=2.0*x(2)+1.0      
    ELSE
      Prim(1)=1.0
      Prim(4)=x(2)+1.5
    END IF


    ro = Prim(1)
    p  = Prim(4)
    ss=SQRT(Gmm*p/ro)
    Prim(2)=0.0
    Prim(3)=-0.025*ss*COS(8.0*PI*x(1))


#ifdef MULTIFLUID
    Prim(5)=0.0
#endif

#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)



  !*------------------------------------------
  !*[13] Shock-vortex interaction
  !*------------------------------------------
  CASE(13) 

    ro1=1.0 
    u1 =SQRT(Gmm)*Ma_shock
    v1 =0.0
    p1 =1.0

    ro2=1.86206896551724133
    u2 =0.95314618727716027
    v2 =0.0                
    p2 =2.45833333333333333

    IF (x(1) .LT. 0.5) THEN
      Prim(1)=ro1
      Prim(2)=u1
      Prim(3)=v1
      Prim(4)=p1
    ELSE 
      Prim(1)=ro2
      Prim(2)=u2
      Prim(3)=v2
      Prim(4)=p2
    END IF

    vm = Ma_vortex*SQRT(Gmm)
    a  = 0.075
    b  = 0.175
    r  = SQRT( (x(1)-0.25)**2 + (x(2)-0.5)**2 )
    T1 = p1/(ro1*R_EOS)

    Ccoeff = T1
    Bcoeff = T1 - (Gmm-1.0)/(R_EOS*Gmm)*vm**2*a**2/(a**2-b**2)**2*(b**2/2.0-2.0*b**2*LOG(b)+b**4*b**(-2)/(-2.0))
    Acoeff = Bcoeff + (Gmm-1.0)/(R_EOS*Gmm)*vm**2*a**2/(a**2-b**2)**2*(a**2/2.0-2.0*b**2*LOG(a)+b**4*a**(-2)/(-2.0))-(Gmm-1.0)/(R_EOS*Gmm)*vm**2/2.0



    IF (r .LE. a) THEN
      vtheta=vm*r/a
      Temp  =Acoeff+(Gmm-1.0)/(R_EOS*Gmm)*vm**2/a**2*r**2/2.0
    ELSEIF (r .GE. b) THEN
      vtheta=0.0
      Temp  =Ccoeff
    ELSE
      vtheta=vm*a/(a**2-b**2)*(r-b**2/r)
      Temp  =Bcoeff+(Gmm-1.0)/(R_EOS*Gmm)*vm**2*a**2/(a**2-b**2)**2*(r**2/2.0-2.0*b**2*LOG(r)+b**4*r**(-2)/(-2.0))
    END IF


    !*Inside vortex only
    IF(r .LT. b) THEN
      x_d(1)=x(1)-0.25
      x_d(2)=x(2)-0.5
      theta = ATAN2(x_d(2), x_d(1))

      Prim(1)=ro1 * (Temp/T1)**(1.0/(Gmm-1.0))
      Prim(2)=Prim(2)-SIN(theta)*vtheta
      Prim(3)=Prim(3)+COS(theta)*vtheta
      Prim(4)=p1  * (Temp/T1)**(Gmm/(Gmm-1.0))

    END IF



#ifdef MULTIFLUID
    Prim(5)=0.0
#endif

#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)




  !*------------------------------------------
  !*[892] Smooth periodic IC with the purpose of verifying conservation
  !*------------------------------------------
  CASE(892)
    Cons(1)=8.0+0.5*SIN(2.0*Pi*x(1))
    Cons(2)=0.5+COS(8.0*Pi*x(2))
    Cons(3)=0.5+SIN(4.0*Pi*x(2))
    Cons(4)=7.0+SIN(4.0*Pi*x(1))

#ifdef MULTIFLUID
    Cons(5)=0.0
#endif

! #ifdef MOMENTUMINPRIMITIVEVARIABLES
!     Prim(2:3)=Prim(2:3)*Prim(1)
! #endif
!     CALL PrimToCons(Prim,Cons)


  !*------------------------------------------
  !*[18] Shock-Tube Problem - 1D in X direction Multifluid
  !*------------------------------------------
  CASE(18) 


    Prim(2)=0.
    Prim(3)=0.

    IF (x(1) .LT. 0.5) THEN
      Prim(1)=1.0
      Prim(4)=1.0
#ifdef MULTIFLUID
      Prim(5)=1.0
#endif
    ELSE
      Prim(1)=0.125
      Prim(4)=0.1
#ifdef MULTIFLUID
      Prim(5)=-1.0
#endif
    END IF

#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)

  !*------------------------------------------
  !*[19] Shock-Tube Problem - 1D in Y direction Multifluid
  !*------------------------------------------
  CASE(19) 


    Prim(2)=0.
    Prim(3)=0.

    IF (x(2) .LE. 0.5) THEN
      Prim(1)=1.0
      Prim(4)=1.0
#ifdef MULTIFLUID
      Prim(5)=1.0
#endif
    ELSE
      Prim(1)=0.125
      Prim(4)=0.1
#ifdef MULTIFLUID
      Prim(5)=-1.0
#endif
    END IF

#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)


  !*------------------------------------------
  !*[28] Shock-Tube Problem - 1D in X direction Multifluid
  !*------------------------------------------
  CASE(28) 


    Prim(2)=0.
    Prim(3)=0.

    IF (x(1) .LE. 0.5) THEN
      Prim(1)=1.0
      Prim(4)=1.0
    ELSE
      Prim(1)=0.125
      Prim(4)=0.1
    END IF

#ifdef MULTIFLUID
      level_set=-(x(1)-0.5)

#ifdef PRIMITIVEFORMULATIONLEVELSET
      Prim(5)=level_set
#else
      Prim(5)=Prim(1)*level_set
#endif

#endif



#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)


  !*------------------------------------------
  !*[-28] Shock-Tube Problem - 1D in X direction Multifluid reverse
  !*------------------------------------------
  CASE(-28) 


    Prim(2)=0.
    Prim(3)=0.

    IF (x(1) .LE. 0.5) THEN
      Prim(1)=0.125
      Prim(4)=0.1  
    ELSE
      Prim(1)=1.0 
      Prim(4)=1.0 
    END IF

#ifdef MULTIFLUID
      level_set=-(x(1)-0.5)

#ifdef PRIMITIVEFORMULATIONLEVELSET
      Prim(5)=level_set
#else
      Prim(5)=Prim(1)*level_set
#endif

#endif



#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)


  !*------------------------------------------
  !*[29] Shock-Tube Problem - 1D in Y direction Multifluid
  !*------------------------------------------
  CASE(29) 


    Prim(2)=0.
    Prim(3)=0.

    IF (x(2) .LE. 0.5) THEN
      Prim(1)=1.0
      Prim(4)=1.0
    ELSE
      Prim(1)=0.125
      Prim(4)=0.1
    END IF

#ifdef MULTIFLUID
      level_set=-(x(2)-0.5)

#ifdef PRIMITIVEFORMULATIONLEVELSET
      Prim(5)=level_set
#else
      Prim(5)=Prim(1)*level_set
#endif

#endif


#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)


  !*------------------------------------------
  !*[30] Stiff Shock-Tube Problem - 1D in X direction Multifluid
  !*------------------------------------------
  CASE(30) 


    Prim(2)=0.
    Prim(3)=0.

    IF (x(1) .LE. 0.5) THEN
      Prim(1)=1.0
      Prim(4)=500.0
    ELSE
      Prim(1)=1.0
      Prim(4)=0.2
    END IF

#ifdef MULTIFLUID
      level_set=-(x(1)-0.5)

#ifdef PRIMITIVEFORMULATIONLEVELSET
      Prim(5)=level_set
#else
      Prim(5)=Prim(1)*level_set
#endif

    ! IF (x(1) .LE. 0.5) THEN
    !   Prim(5)=1.0*Prim(1)
    ! ELSE
    !   Prim(5)=-1.0*Prim(1)
    ! END IF

#endif



#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)

  !*------------------------------------------
  !*[31] Water-Air Model Using the Stiff Equation of State Multifluid
  !*------------------------------------------
  CASE(31) 

    Prim(2)=0.
    Prim(3)=0.

    IF (x(1) .LE. 0.7) THEN
      Prim(1)=1000.0
      Prim(4)=10.0**9
    ELSE
      Prim(1)=50.0
      Prim(4)=10.0**5
    END IF

#ifdef MULTIFLUID
      level_set=-(x(1)-0.7)

#ifdef PRIMITIVEFORMULATIONLEVELSET
      Prim(5)=level_set
#else
      Prim(5)=Prim(1)*level_set
#endif

#endif



#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)

  !*------------------------------------------
  !*[40] Helium Bubble Multifluid
  !*------------------------------------------
  CASE(40)     

    r=SQRT(x(1)**2+x(2)**2)

    Prim(3)=0.

    IF (x(1) .GT. 0.75) THEN !*ZONE C
      Prim(1)=4.0/3.0
      Prim(2)=-0.3535
      Prim(4)=1.5
#ifdef MULTIFLUID
      level_set=-1.0
      Prim(5)=level_set*Prim(1)
#endif
    ELSE IF (r .LT. 0.25) THEN  !*ZONE A
      Prim(1)=4.0/29.0
      Prim(2)=0.
      Prim(4)=1.0
#ifdef MULTIFLUID
      level_set=1.0
      Prim(5)=level_set*Prim(1)
#endif
    ELSE  !*ZONE B
      Prim(1)=1.0
      Prim(2)=0.
      Prim(4)=1.0
#ifdef MULTIFLUID
      level_set=-1.0
      Prim(5)=level_set*Prim(1)
#endif
    END IF


#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)

  !*------------------------------------------
  !*[41] R22 Bubble Multifluid
  !*------------------------------------------
  CASE(41)     

    r=SQRT(x(1)**2+x(2)**2)

    Prim(3)=0.

    IF (x(1) .GT. 0.75) THEN !*ZONE C
      Prim(1)=4.0/3.0
      Prim(2)=-0.3535
      Prim(4)=1.5
#ifdef MULTIFLUID
      level_set=-1.0
      Prim(5)=level_set*Prim(1)
#endif
    ELSE IF (r .LT. 0.25) THEN  !*ZONE A
      Prim(1)=3.1538
      Prim(2)=0.
      Prim(4)=1.0
#ifdef MULTIFLUID
      level_set=1.0
      Prim(5)=level_set*Prim(1)
#endif
    ELSE  !*ZONE B
      Prim(1)=1.0
      Prim(2)=0.
      Prim(4)=1.0
#ifdef MULTIFLUID
      level_set=-1.0
      Prim(5)=level_set*Prim(1)
#endif
    END IF


#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)

  CASE(-300) ! Shock-turbulence interaction Shu-Osher
    !*Efficient implementation of essentially non-oscillatory shock capturing schemes, II, Shu, Osher
    IF(x(1) .LT. -4.0) THEN
      Prim(1) = 3.857143      
      Prim(2) = 2.629369
      Prim(3) = 0.0
      Prim(4) = 10.333333
    ELSE
      Prim(1) = 1.0+0.2*SIN(5.0*x(1))     
      Prim(2) = 0.      
      Prim(3) = 0.0
      Prim(4) = 1.0      
    END IF

#ifdef MULTIFLUID
      level_set=0.0
      Prim(5)=level_set*Prim(1)
#endif


#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)

  CASE(-301) ! Shock-turbulence interaction Shu-Osher modified by Titarev and Toro
    !*Finite-volume WENO schemes for three-dimensional conservation laws, V.A. Titarev, E.F. Toro 
    IF(x(1) .LT. -4.5) THEN
      Prim(1) = 1.515695      
      Prim(2) = 0.523346
      Prim(3) = 0.0
      Prim(4) = 1.80500
    ELSE
      Prim(1) = 1.0+0.1*SIN(20.0*PI*x(1))     
      Prim(2) = 0.0      
      Prim(3) = 0.0
      Prim(4) = 1.0      
    END IF

#ifdef MULTIFLUID
      level_set=0.0
      Prim(5)=level_set*Prim(1)
#endif


#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)


  CASE(-302,-304) ! Shock-turbulence interaction Shu-Osher modified example 1
    IF(x(1) .LT. -4.5) THEN
      Prim(1) = 1.515695      
      Prim(2) = 0.523346
      Prim(3) = 0.0
      Prim(4) = 1.80500
    ELSE
      Prim(1) = 1.0+0.1*SIN(20.0*x(1))     
      Prim(2) = 0.0      
      Prim(3) = 0.0
      Prim(4) = 1.0      
    END IF

#ifdef MULTIFLUID
      level_set=0.0
      Prim(5)=level_set*Prim(1)
#endif


#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)


  CASE(-303,-305) ! Shock-turbulence interaction Shu-Osher modified example 2
    IF(x(1) .LT. -4.0) THEN
      Prim(1) = 27.0/7.0      
      Prim(2) = 4.0/9.0*SQRT(35.0)
      Prim(3) = 0.0
      Prim(4) = 31.0/3.0
    ELSE
      Prim(1) = 1.0+0.2*SIN(5.0*x(1))     
      Prim(2) = 0.0      
      Prim(3) = 0.0
      Prim(4) = 1.0      
    END IF

#ifdef MULTIFLUID
      level_set=0.0
      Prim(5)=level_set*Prim(1)
#endif


#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)


  CASE(350) ! Woodword-Colella 
    Prim(1) = 1.0
    Prim(2) = 0.0
    IF(x(1).LT. 0.1) THEN
      Prim(4) = 1000.0    
    ELSE IF (x(1).GT. 0.9) THEN
      Prim(4) = 100.0    
    ELSE
      Prim(4) = 1e-2    
    END IF

    Prim(3)=0.0
#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)

  CASE(401) ! RP1 
    IF(x(1).LT. DiscontinuityLocation) THEN
      Prim(1) = 1.0
      Prim(2) = 0.75
      Prim(4) = 1.0
    ELSE
      Prim(1) = 0.125
      Prim(2) = 0.0
      Prim(4) = 0.1
    END IF

    Prim(3)=0.0
#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)

  CASE(402) ! RP2
    IF(x(1).LT. DiscontinuityLocation) THEN
      Prim(1) = 1.0
      Prim(2) = -2.0
      Prim(4) = 0.4
    ELSE
      Prim(1) = 1.0
      Prim(2) = 2.0
      Prim(4) = 0.4
    END IF

    Prim(3)=0.0
#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)

  CASE(412) ! RP2 !*Modified velocity
    IF(x(1).LT. DiscontinuityLocation) THEN
      Prim(1) = 1.0
      Prim(2) = -1.0
      Prim(4) = 0.4
    ELSE
      Prim(1) = 1.0
      Prim(2) = 1.0
      Prim(4) = 0.4
    END IF

    Prim(3)=0.0
#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)

  CASE(403) ! RP3
    IF(x(1).LT. DiscontinuityLocation) THEN
      Prim(1) = 1.0
      Prim(2) = 0.0
      Prim(4) = 1000.0
    ELSE
      Prim(1) = 1.0
      Prim(2) = 0.0
      Prim(4) = 0.01
    END IF

    Prim(3)=0.0
#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)

  CASE(404) ! RP4
    IF(x(1).LT. DiscontinuityLocation) THEN
      Prim(1) = 5.99924
      Prim(2) = 19.5975
      Prim(4) = 460.894
    ELSE
      Prim(1) = 5.99242
      Prim(2) = -6.19633
      Prim(4) = 46.0950
    END IF

    Prim(3)=0.0
#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)

  CASE(405) ! RP5
    IF(x(1).LT. DiscontinuityLocation) THEN
      Prim(1) = 1.0
      Prim(2) = -19.59745
      Prim(4) = 1000.0
    ELSE
      Prim(1) = 1.0
      Prim(2) = -19.59745
      Prim(4) = 0.01
    END IF

    Prim(3)=0.0
#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)

  CASE(406) ! RP6
    IF(x(1).LT. DiscontinuityLocation) THEN
      Prim(1) = 1.0
      Prim(2) = 2.0
      Prim(4) = 0.1
    ELSE
      Prim(1) = 1.0
      Prim(2) = -2.0
      Prim(4) = 0.1
    END IF

    Prim(3)=0.0
#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)

  CASE(420) ! Stationary contact
    IF(x(1) .LT. DiscontinuityLocation) THEN
      Prim(1) = 1.4
      Prim(2) = 0.0
      Prim(3) = 0.0
      Prim(4) = 1.0
    ELSE
      Prim(1) = 1.0
      Prim(2) = 0.0
      Prim(3) = 0.0
      Prim(4) = 1.0
    END IF

#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)

  CASE(421) ! Moving contact
    IF(x(1) .LT. DiscontinuityLocation) THEN
      Prim(1) = 1.4
      Prim(2) = 0.1
      Prim(3) = 0.0
      Prim(4) = 1.0
    ELSE
      Prim(1) = 1.0
      Prim(2) = 0.1
      Prim(3) = 0.0
      Prim(4) = 1.0
    END IF

#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)

  CASE(500,501) ! Lax
    IF (x(1).LE. 5.0) THEN
      Prim(1)=0.445*Lambda_Lax 
      Prim(2)=0.698
      Prim(4)=3.528*Lambda_Lax
    ELSE
      Prim(1)=0.5*Lambda_Lax 
      Prim(2)=0.0
      Prim(4)=0.571*Lambda_Lax
    END IF


    Prim(3)=0.0
#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)





#ifdef RELAXATION
  !*------------------------------------------
  !*[99,100] Gresho vortex longer time and Gresho vortex
  !*------------------------------------------
  CASE(-99,99,100) 

    !*Center of the vortex
    xc=0.5*(MESH_X1+MESH_X0)

    !*Distance squared from the center of the vortex
    x_d=x-xc
    r2=x_d(1)**2+x_d(2)**2
    r =SQRT(r2)

    !*Angular position    
    theta=ATAN2(x_d(2),x_d(1))

    sintheta=SIN(theta)
    costheta=COS(theta)

    !*Parameters
    ro_bar  = 1.0
    vmax    = 1.0
    p0      = 0.5

    !*Parameters
    Prim(1) = ro_bar

    IF ( r .LT. 0.2 ) THEN
      Omega   = 5.0*r
      Prim(4) = p0 + (25.0/2.0*r2)*EPS_LM**2
    ELSE IF ( (r .GE. 0.2) .AND. (r .LT. 0.4) ) THEN
      Omega   = 2.0 - 5.0*r
      Prim(4) = p0 + (25.0/2.0*r2 + 4.0 * (1.0 - 5.0*r- LOG(0.2) + LOG(r)))*EPS_LM**2
    ELSE
      Omega   = 0.0
      Prim(4) = p0 +(- 2.0 + 4.0*LOG(2.0))*EPS_LM**2
    END IF

    Prim(2)=omega*(-sintheta)
    Prim(3)=omega*( costheta)

#ifdef MULTIFLUID
    Prim(5)=0.0
#endif

#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)

  !*------------------------------------------
  !*[101] Smooth Gresho vortex
  !*------------------------------------------
  CASE(101) 

    !*Center of the vortex
    xc=0.5*(MESH_X1+MESH_X0)

    !*Distance squared from the center of the vortex
    x_d=x-xc
    r2=x_d(1)**2+x_d(2)**2
    r =SQRT(r2)


    !*Angular position    
    theta=ATAN2(x_d(2),x_d(1))

    sintheta=SIN(theta)
    costheta=COS(theta)


    !*Parameters
    ro_bar  = 1.0
    vmax    = 1.0
    p0      = 0.5

    !*Parameters
    Prim(1) = ro_bar

    IF ( r .LT. 0.2 ) THEN
      Omega   = 75.0*r2-250.0*r**3
      Prim(4) = p0 + (1406.25*r**4-7500.0*r**5+(10416.0+2.0/3.0)*r**6)*EPS_LM**2
    ELSE IF ( (r .GE. 0.2) .AND. (r .LT. 0.4) ) THEN
      Omega   = -4.0+60.0*r-225.0*r2+250.0*r**3
      p2=pfunction(r)
      Prim(4) = p0 + (p2)*EPS_LM**2
    ELSE
      Omega   = 0.0
      p2=pfunction(0.4)
      Prim(4) = p0 + (p2)*EPS_LM**2
    END IF

    Prim(2)=omega*(-sintheta)
    Prim(3)=omega*( costheta)

#ifdef MULTIFLUID
    Prim(5)=0.0
#endif

#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)


  !*------------------------------------------
  !*[200] 2d Low-Mach Vortex
  !*------------------------------------------
  CASE(200) 

    u_inf=1.0
    v_inf=1.0

    !*Original center of the vortex before the moevement
    xc=0.5*(MESH_X1+MESH_X0)

    !*We want to get the initial position of x before the movement

    !*x with respect to bottom-left corner
    x_wrt_BL=x-MESH_X0

    !*x with respect to bottom-left corner before movement
    x_wrt_BL_bm(1)=x_wrt_BL(1)-u_inf*t
    x_wrt_BL_bm(2)=x_wrt_BL(2)-v_inf*t

    !*This is the position before movement modulo the length of the domain
    !*NB: MODULO RESULT IS ALWAYS POSITIVE
    x_wrt_BL_bm(1)=MODULO( x_wrt_BL_bm(1), MESH_SX(1) )
    x_wrt_BL_bm(2)=MODULO( x_wrt_BL_bm(2), MESH_SX(2) )

    !*This is the initial position
    x_0=MESH_X0+x_wrt_BL_bm

    !*Distance squared from the center of the vortex at the initial time
    x_d=x_0-xc

    r2=x_d(1)**2+x_d(2)**2

    delta_ro=-(Gmm-1.0)/(8.0*Gmm*Pi**2)*EPS_LM**2*EXP(1.0-r2)
    delta_u=EPS_LM*EXP(0.5*(1.0-r2))/(2.0*Pi)*( -x_d(2) )
    delta_v=EPS_LM*EXP(0.5*(1.0-r2))/(2.0*Pi)*x_d(1)

    Prim(1)=(1.0+delta_ro)**( 1.0 / (Gmm-1.0) )
    Prim(2)=u_inf+delta_u
    Prim(3)=v_inf+delta_v

    delta_energy=EPS_LM**2*( Prim(1)**Gmm/(Gmm-1.0)+0.5*Prim(1)*(Prim(2)**2+Prim(3)**2) )

    Cons(4)=1.0+delta_energy
    Cons(1)=Prim(1)
    Cons(2)=Prim(1)*Prim(2)
    Cons(3)=Prim(1)*Prim(3)

#ifdef MULTIFLUID
    Prim(5)=0.0
    Cons(5)=0.0
#endif

  !*------------------------------------------
  !*[201] 2d Low-Mach Vortex modified
  !*------------------------------------------
  CASE(201) 

    u_inf=20.0*EPS_LM
    v_inf=20.0*EPS_LM

    !*Original center of the vortex before the moevement
    xc=0.5*(MESH_X1+MESH_X0)

    !*We want to get the initial position of x before the movement

    !*x with respect to bottom-left corner
    x_wrt_BL=x-MESH_X0

    !*x with respect to bottom-left corner before movement
    x_wrt_BL_bm(1)=x_wrt_BL(1)-u_inf*t
    x_wrt_BL_bm(2)=x_wrt_BL(2)-v_inf*t

    !*This is the position before movement modulo the length of the domain
    !*NB: MODULO RESULT IS ALWAYS POSITIVE
    x_wrt_BL_bm(1)=MODULO( x_wrt_BL_bm(1), MESH_SX(1) )
    x_wrt_BL_bm(2)=MODULO( x_wrt_BL_bm(2), MESH_SX(2) )

    !*This is the initial position
    x_0=MESH_X0+x_wrt_BL_bm

    !*Distance squared from the center of the vortex at the initial time
    x_d=x_0-xc

    r2=x_d(1)**2+x_d(2)**2

    delta_ro=-(Gmm-1.0)/(8.0*Gmm*Pi**2)*EPS_LM**2*EXP(1.0-r2)
    delta_u=EPS_LM*EXP(0.5*(1.0-r2))/(2.0*Pi)*( -x_d(2) )
    delta_v=EPS_LM*EXP(0.5*(1.0-r2))/(2.0*Pi)*x_d(1)

    Prim(1)=(1.0+delta_ro)**( 1.0 / (Gmm-1.0) )
    Prim(2)=u_inf+delta_u
    Prim(3)=v_inf+delta_v

    delta_energy=EPS_LM**2*( Prim(1)**Gmm/(Gmm-1.0)+0.5*Prim(1)*(Prim(2)**2+Prim(3)**2) )

    Cons(4)=1.0+delta_energy
    Cons(1)=Prim(1)
    Cons(2)=Prim(1)*Prim(2)
    Cons(3)=Prim(1)*Prim(3)

#ifdef MULTIFLUID
    Prim(5)=0.0
    Cons(5)=0.0
#endif

  !*------------------------------------------
  !*[300] Double Shear Layer
  !*------------------------------------------
  !*A novel full-Euler low Mach number IMEX splitting, Jonas Zeifang, Jochen Schütz, Klaus Kaiser, Andrea Beck,Maria Lukacova-Medvidova, Sebastian Noelle3
  !*------------------------------------------
  CASE(300) 


    Prim(1)=PI/15.0
    IF (x(2) .LE. PI) THEN
      Prim(2)=TANH( ( x(2)-0.5*PI )/Prim(1) )
    ELSE
      Prim(2)=TANH( ( 1.5*PI-x(2) )/Prim(1) )
    END IF
    Prim(3)=0.05*SIN(x(1))
    Prim(4)=1.0/Gmm


#ifdef MULTIFLUID
    Prim(5)=0.0
#endif

#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)

  !*------------------------------------------
  !*[301] Baroclinic vorticity generation problem
  !*------------------------------------------
  !*A novel full-Euler low Mach number IMEX splitting, Jonas Zeifang, Jochen Schütz, Klaus Kaiser, Andrea Beck,Maria Lukacova-Medvidova, Sebastian Noelle3
  !*------------------------------------------
  CASE(301) 


    ro0 = 1.0
    u0  = SQRT(Gmm)
    p0  = 1.0

    Prim(1)=ro0+EPS_LM/(2000.0)*(1.0+COS( PI*x(1)*EPS_LM )) + phifunction(x(2))
    Prim(2)=0.5*u0*(1.0+COS( PI*x(1)*EPS_LM ))
    Prim(3)=0.0
    Prim(4)=p0+0.5*EPS_LM*Gmm*(1.0+COS( PI*x(1)*EPS_LM ))

#ifdef MULTIFLUID
    Prim(5)=0.0
#endif

#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)


  !*------------------------------------------
  !*[302] Kelvin-Helmholtz instability
  !*------------------------------------------
  !*All-speed numerical methods for the Euler equations via a sequential explicit time integration, Wasilij Barsukow
  !*------------------------------------------
  CASE(302) 

    IF (x(2) .LE. 0.5) THEN
      Prim(1)=1.001
      Prim(2)=+0.1
    ELSE
      Prim(1)=0.999
      Prim(2)=-0.1
    END IF
    Prim(3)=0.001*SIN(2.0*PI*x(1))
    Prim(4)=5.0

#ifdef MULTIFLUID
    Prim(5)=0.0
#endif

#ifdef MOMENTUMINPRIMITIVEVARIABLES
    Prim(2:3)=Prim(2:3)*Prim(1)
#endif
    CALL PrimToCons(Prim,Cons)


#endif



#endif
  CASE DEFAULT
    ErrorMessage = "Exact function not specified"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
CONTAINS
   
#ifdef SW
  REAL FUNCTION hSmoothAuxiliary(x)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x

    hSmoothAuxiliary=1.-0.5*exp(-1./atan(1.-x)**3.)

  END FUNCTION

  REAL FUNCTION hDerivSmoothAuxiliary(x)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x

    hDerivSmoothAuxiliary=3.*0.5*exp(1./atan(x - 1.)**3.)/(atan(x - 1.)**4*((x - 1.)**2 + 1.))

  END FUNCTION
#endif

  REAL FUNCTION pfunction(r)
    IMPLICIT NONE
    REAL, INTENT(IN) :: r

    pfunction=65.8843399322788 - 480.0*r + 2700.0*r**2 - (9666.0+2.0/3.0)*r**3 &
             &+20156.25*r**4 - 22500.0*r**5 + (10416.0 + 2.0/3.0)*r**6 + 16.0*LOG(r)

  END FUNCTION

#ifdef RELAXATION
  REAL FUNCTION phifunction(y)
    IMPLICIT NONE
    REAL, INTENT(IN) :: y
    REAL             :: Ly

    Ly=2.0/5.0*(1.0/EPS_LM)

    IF( (y .GE. 0.0) .AND. (y .LE. 0.5*Ly) ) THEN
      phifunction=1.8*y/Ly
    ELSE
      phifunction=1.8*(y/Ly-1.0)
    END IF

  END FUNCTION
#endif

END SUBROUTINE ExactFunction
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ExactFunctionWB(WhichInitialCondition,x,Cons)
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: PI 
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DENSITY
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
INTEGER,INTENT(IN) :: WhichInitialCondition
REAL,INTENT(IN)    :: x(1:nDims)
REAL,INTENT(OUT)   :: Cons(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: Prim(1:nVar)
CHARACTER(LEN=255) :: ErrorMessage

SELECT CASE (WhichInitialCondition)

  CASE DEFAULT
    ErrorMessage = "Exact WB function not specified"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

END SUBROUTINE ExactFunctionWB
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE SourceTerms(t)
!-------------------------------------------------------------------------------!
USE MOD_Reconstruction     ,ONLY: MUSCL
USE MOD_Reconstruction     ,ONLY: WENO3_SecondSweep 
USE MOD_Reconstruction     ,ONLY: WENO5_SecondSweep 
USE MOD_Reconstruction     ,ONLY: WENO7_SecondSweep 
USE MOD_FiniteVolume2D_vars,ONLY: S
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP  
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGP
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: GravitationalPotentialFlag
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: S_in_qp(1:nVar,nGPs,nGPs,nElemsX,nElemsY) !*Source in quadrature points
REAL             :: Vtemp(1:nVar,-nGhosts:nElemsX+nGhosts+1,1:nElemsY,1:nGPs)
REAL             :: Vtemp2(1:nVar,1:nGPs,1:nGPs,1:nElemsX,1:nElemsY)
INTEGER          :: ii, jj, iVar, iGP, jGP
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

S = 0.0 ! S(1:nVar,nElemsX,nElemsY)

!int_{iixjj} S(1:nVar,ii,jj) dxdy

IF (GravitationalPotentialFlag .GT. 0) THEN


  SELECT CASE (ABS(Reconstruction))
    CASE(1,10,2,20,21,22,23,24,25,26,27,28,29)
      DO jj=1,nElemsY
        DO ii=1,nElemsX
          S_in_qp(1:nVar,nGPs,nGPs,ii,jj) = SourceFunc( U(1:nVar,ii,jj) , MeshBary(:,ii,jj) )
        END DO
      END DO
    CASE(3)
      PRINT*, "No option for reconstructing a specific kind of variable in source"
      STOP
      DO iVar=1,nVar
        DO jj=1,nElemsY
          DO ii=-nGhosts,nElemsX+nGhosts+1
             CALL WENO3_SecondSweep( U(iVar,ii,jj-nGhosts:jj+nGhosts) , Vtemp(iVar,ii,jj,1:nGPs) )
          END DO
        END DO
      END DO
      DO jj=1,nElemsY
        DO ii=1,nElemsX
          DO iGP=1,nGPs
            DO iVar=1,nVar
              CALL WENO3_SecondSweep( Vtemp(iVar,ii-nGhosts:ii+nGhosts,jj,iGP) , Vtemp2(iVar,1:nGPs,iGP,ii,jj) )
            END DO
            DO jGP=1,nGPs
              S_in_qp(1:nVar,jGP,iGP,ii,jj) = SourceFunc( Vtemp2(1:nVar,jGP,iGP,ii,jj) , MeshGP(:,ii,jj,jGP,iGP)  )
            END DO 
          END DO 
        END DO
      END DO
    CASE(4,5)
      PRINT*, "No option for reconstructing a specific kind of variable in source"
      STOP
      DO iVar=1,nVar
        DO jj=1,nElemsY
          DO ii=-nGhosts,nElemsX+nGhosts+1
             CALL WENO5_SecondSweep( U(iVar,ii,jj-nGhosts:jj+nGhosts) , Vtemp(iVar,ii,jj,1:nGPs) )
          END DO
        END DO
      END DO
      DO jj=1,nElemsY
        DO ii=1,nElemsX
          DO iGP=1,nGPs
            DO iVar=1,nVar
              CALL WENO5_SecondSweep( Vtemp(iVar,ii-nGhosts:ii+nGhosts,jj,iGP) , Vtemp2(iVar,1:nGPs,iGP,ii,jj) )
            END DO
            DO jGP=1,nGPs
              S_in_qp(1:nVar,jGP,iGP,ii,jj) = SourceFunc( Vtemp2(1:nVar,jGP,iGP,ii,jj) , MeshGP(:,ii,jj,jGP,iGP)  )
            END DO 
          END DO 
        END DO
      END DO
    CASE(7)
      PRINT*, "No option for reconstructing a specific kind of variable in source"
      STOP
      DO iVar=1,nVar
        DO jj=1,nElemsY
          DO ii=-nGhosts,nElemsX+nGhosts+1
             CALL WENO7_SecondSweep( U(iVar,ii,jj-nGhosts:jj+nGhosts) , Vtemp(iVar,ii,jj,1:nGPs) )
          END DO
        END DO
      END DO
      DO jj=1,nElemsY
        DO ii=1,nElemsX
          DO iGP=1,nGPs
            DO iVar=1,nVar
              CALL WENO7_SecondSweep( Vtemp(iVar,ii-nGhosts:ii+nGhosts,jj,iGP) , Vtemp2(iVar,1:nGPs,iGP,ii,jj) )
            END DO
            DO jGP=1,nGPs
              S_in_qp(1:nVar,jGP,iGP,ii,jj) = SourceFunc( Vtemp2(1:nVar,jGP,iGP,ii,jj) , MeshGP(:,ii,jj,jGP,iGP)  )
            END DO 
          END DO 
        END DO
      END DO

    CASE DEFAULT
      ErrorMessage = "Reconstruction not implemented"
      WRITE(*,*) ErrorMessage
      STOP
  END SELECT


  DO jj=1,nElemsY
    DO ii=1,nElemsX
      DO iGP=1,nGPs
        DO jGP=1,nGPs
          S(1:nVar,ii,jj) = S(1:nVar,ii,jj) + WeightsGP(iGP,jGP) * S_in_qp(1:nVar,iGP,jGP,ii,jj) 
        END DO 
      END DO 
    END DO
  END DO

END IF
!-------------------------------------------------------------------------------!
END SUBROUTINE SourceTerms
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE SourceTerm_Conserved_INPUT_CONSERVED(t,W,vbXs,vbXe,vbYs,vbYe,lbXs,lbXe,lbYs,lbYe,MB,MBbXs,MBbXe,MBbYs,MBbYe,MGP,MGPbXs,MGPbXe,MGPbYs,MGPbYe)
!-------------------------------------------------------------------------------!
USE MOD_Reconstruction     ,ONLY: MUSCL
USE MOD_Reconstruction     ,ONLY: WENO3_SecondSweep 
USE MOD_Reconstruction     ,ONLY: WENO5_SecondSweep 
USE MOD_Reconstruction     ,ONLY: WENO7_SecondSweep 
USE MOD_FiniteVolume2D_vars,ONLY: S
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGP
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: GravitationalPotentialFlag
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: t
REAL,                                        INTENT(IN)           :: W(1:nVar,vbXs:vbXe,vbYs:vbYe)
INTEGER,                                     INTENT(IN)           :: vbXs,vbXe,vbYs,vbYe     !*Bounds for the vector
INTEGER,                                     INTENT(IN)           :: lbXs,lbXe,lbYs,lbYe     !*Bounds for the loop   
REAL   ,                                     INTENT(IN)           :: MB(1:nDims,MBbXs:MBbXe,MBbYs:MBbYe)
INTEGER,                                     INTENT(IN)           :: MBbXs,MBbXe,MBbYs,MBbYe !*Bounds for mesh bary for the loop   
REAL   ,                                     INTENT(IN)           :: MGP(1:nDims,MGPbXs:MGPbXe,MGPbYs:MGPbYe,1:nGPs,1:nGPs)
INTEGER,                                     INTENT(IN)           :: MGPbXs,MGPbXe,MGPbYs,MGPbYe !*Bounds for mesh bary for the loop   
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: S_in_qp(1:nVar,nGPs,nGPs,nElemsX,nElemsY) !*Source in quadrature points
REAL             :: Vtemp(1:nVar,-nGhosts:nElemsX+nGhosts+1,1:nElemsY,1:nGPs)
REAL             :: Vtemp2(1:nVar,1:nGPs,1:nGPs,1:nElemsX,1:nElemsY)
INTEGER          :: ii, jj, iVar, iGP, jGP
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

S = 0.0 ! S(1:nVar,nElemsX,nElemsY)

!int_{iixjj} S(1:nVar,ii,jj) dxdy

IF (GravitationalPotentialFlag .GT. 0) THEN

  SELECT CASE (ABS(Reconstruction))
    CASE(1,10,2,20,21,22,23,24,25,26,27,28,29)
      DO jj=1,nElemsY
        DO ii=1,nElemsX
          S_in_qp(1:nVar,nGPs,nGPs,ii,jj) = SourceFunc( W(1:nVar,ii,jj) , MB(:,ii,jj) )
        END DO
      END DO
    CASE(3)
      PRINT*, "No option for reconstructing a specific kind of variable in source"
      STOP
      DO iVar=1,nVar
        DO jj=1,nElemsY
          DO ii=-nGhosts,nElemsX+nGhosts+1
             CALL WENO3_SecondSweep( W(iVar,ii,jj-nGhosts:jj+nGhosts) , Vtemp(iVar,ii,jj,1:nGPs) )
          END DO
        END DO
      END DO
      DO jj=1,nElemsY
        DO ii=1,nElemsX
          DO iGP=1,nGPs
            DO iVar=1,nVar
              CALL WENO3_SecondSweep( Vtemp(iVar,ii-nGhosts:ii+nGhosts,jj,iGP) , Vtemp2(iVar,1:nGPs,iGP,ii,jj) )
            END DO
            DO jGP=1,nGPs
              S_in_qp(1:nVar,jGP,iGP,ii,jj) = SourceFunc( Vtemp2(1:nVar,jGP,iGP,ii,jj) , MGP(:,ii,jj,jGP,iGP)  )
            END DO 
          END DO 
        END DO
      END DO
    CASE(4,5)
      PRINT*, "No option for reconstructing a specific kind of variable in source"
      STOP
      DO iVar=1,nVar
        DO jj=1,nElemsY
          DO ii=-nGhosts,nElemsX+nGhosts+1
             CALL WENO5_SecondSweep( W(iVar,ii,jj-nGhosts:jj+nGhosts) , Vtemp(iVar,ii,jj,1:nGPs) )
          END DO
        END DO
      END DO
      DO jj=1,nElemsY
        DO ii=1,nElemsX
          DO iGP=1,nGPs
            DO iVar=1,nVar
              CALL WENO5_SecondSweep( Vtemp(iVar,ii-nGhosts:ii+nGhosts,jj,iGP) , Vtemp2(iVar,1:nGPs,iGP,ii,jj) )
            END DO
            DO jGP=1,nGPs
              S_in_qp(1:nVar,jGP,iGP,ii,jj) = SourceFunc( Vtemp2(1:nVar,jGP,iGP,ii,jj) , MGP(:,ii,jj,jGP,iGP)  )
            END DO 
          END DO 
        END DO
      END DO
    CASE(7)
      PRINT*, "No option for reconstructing a specific kind of variable in source"
      STOP
      DO iVar=1,nVar
        DO jj=1,nElemsY
          DO ii=-nGhosts,nElemsX+nGhosts+1
             CALL WENO7_SecondSweep( W(iVar,ii,jj-nGhosts:jj+nGhosts) , Vtemp(iVar,ii,jj,1:nGPs) )
          END DO
        END DO
      END DO
      DO jj=1,nElemsY
        DO ii=1,nElemsX
          DO iGP=1,nGPs
            DO iVar=1,nVar
              CALL WENO7_SecondSweep( Vtemp(iVar,ii-nGhosts:ii+nGhosts,jj,iGP) , Vtemp2(iVar,1:nGPs,iGP,ii,jj) )
            END DO
            DO jGP=1,nGPs
              S_in_qp(1:nVar,jGP,iGP,ii,jj) = SourceFunc( Vtemp2(1:nVar,jGP,iGP,ii,jj) , MGP(:,ii,jj,jGP,iGP)  )
            END DO 
          END DO 
        END DO
      END DO

    CASE DEFAULT
      ErrorMessage = "Reconstruction not implemented"
      WRITE(*,*) ErrorMessage
      STOP
  END SELECT


  DO jj=1,nElemsY
    DO ii=1,nElemsX
      DO iGP=1,nGPs
        DO jGP=1,nGPs
          S(1:nVar,ii,jj) = S(1:nVar,ii,jj) + WeightsGP(iGP,jGP) * S_in_qp(1:nVar,iGP,jGP,ii,jj) 
        END DO 
      END DO 
    END DO
  END DO

END IF
!-------------------------------------------------------------------------------!
END SUBROUTINE SourceTerm_Conserved_INPUT_CONSERVED
!===============================================================================!
!
!
!
!===============================================================================!
FUNCTION SourceFunc(Q,X) RESULT(S) !*ALERT, it should be in conserved variables but I never tested it
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nVar
IMPLICIT NONE
REAL, DIMENSION(1:nVar), INTENT(IN)  :: Q 
REAL, DIMENSION(1:nDims) , INTENT(IN)  :: X
REAL, DIMENSION(1:nVar) :: S 

S(1) = 0.
S(2) = -Q(1)*Gravitational_Potential_X(X) 
S(3) = -Q(1)*Gravitational_Potential_Y(X) 
S(4) = -(Q(2)*Gravitational_Potential_X(X)+Q(3)*Gravitational_Potential_Y(X))

!-------------------------------------------------------------------------------!
END FUNCTION SourceFunc
!===============================================================================!
!
!
!
!===============================================================================!
REAL FUNCTION Gravitational_Potential(X)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: PI   
USE MOD_FiniteVolume2D_vars,ONLY: GravitationalPotentialFlag   
IMPLICIT NONE
REAL, DIMENSION(1:nDims) , INTENT(IN)  :: X
REAL                            :: r2

SELECT CASE (GravitationalPotentialFlag)
  CASE(1)
    Gravitational_Potential = -X(2)
  CASE DEFAULT
    Gravitational_Potential = 0.
END SELECT

!-------------------------------------------------------------------------------!
END FUNCTION Gravitational_Potential
!===============================================================================!
!
!
!
!===============================================================================!
REAL FUNCTION Gravitational_Potential_X(X)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: PI   
USE MOD_FiniteVolume2D_vars,ONLY: GravitationalPotentialFlag   
IMPLICIT NONE
REAL, DIMENSION(1:nDims) , INTENT(IN)  :: X
REAL                            :: r2

SELECT CASE (GravitationalPotentialFlag)
  CASE(1)
    Gravitational_Potential_X = 0.
  CASE DEFAULT
    Gravitational_Potential_X = 0.
END SELECT

!-------------------------------------------------------------------------------!
END FUNCTION Gravitational_Potential_X
!===============================================================================!
!
!
!
!===============================================================================!
REAL FUNCTION Gravitational_Potential_Y(X)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: PI   
USE MOD_FiniteVolume2D_vars,ONLY: GravitationalPotentialFlag   
IMPLICIT NONE
REAL, DIMENSION(1:nDims) , INTENT(IN)  :: X
REAL                            :: r2

SELECT CASE (GravitationalPotentialFlag)
  CASE(1)
    Gravitational_Potential_Y = -1.0
  CASE DEFAULT
    Gravitational_Potential_Y = 0.
END SELECT

!-------------------------------------------------------------------------------!
END FUNCTION Gravitational_Potential_Y
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE SourceTerm_Primitive_INPUT_PRIMITIVE(t,W,vbXs,vbXe,vbYs,vbYe,lbXs,lbXe,lbYs,lbYe,MB,MBbXs,MBbXe,MBbYs,MBbYe,MGP,MGPbXs,MGPbXe,MGPbYs,MGPbYe)
!-------------------------------------------------------------------------------!
USE MOD_Reconstruction     ,ONLY: MUSCL
USE MOD_Reconstruction     ,ONLY: WENO3_SecondSweep 
USE MOD_Reconstruction     ,ONLY: WENO5_SecondSweep 
USE MOD_Reconstruction     ,ONLY: WENO7_SecondSweep 
USE MOD_FiniteVolume2D_vars,ONLY: S
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGP
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: GravitationalPotentialFlag
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: t
REAL,                                        INTENT(IN)           :: W(1:nVar,vbXs:vbXe,vbYs:vbYe)
INTEGER,                                     INTENT(IN)           :: vbXs,vbXe,vbYs,vbYe     !*Bounds for the vector
INTEGER,                                     INTENT(IN)           :: lbXs,lbXe,lbYs,lbYe     !*Bounds for the loop   
REAL   ,                                     INTENT(IN)           :: MB(1:nDims,MBbXs:MBbXe,MBbYs:MBbYe)
INTEGER,                                     INTENT(IN)           :: MBbXs,MBbXe,MBbYs,MBbYe !*Bounds for mesh bary for the loop   
REAL   ,                                     INTENT(IN)           :: MGP(1:nDims,MGPbXs:MGPbXe,MGPbYs:MGPbYe,1:nGPs,1:nGPs)
INTEGER,                                     INTENT(IN)           :: MGPbXs,MGPbXe,MGPbYs,MGPbYe !*Bounds for mesh bary for the loop   
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: S_in_qp(1:nVar,nGPs,nGPs,nElemsX,nElemsY) !*Source in quadrature points
REAL             :: Vtemp(1:nVar,-nGhosts:nElemsX+nGhosts+1,1:nElemsY,1:nGPs)
REAL             :: Vtemp2(1:nVar,1:nGPs,1:nGPs,1:nElemsX,1:nElemsY)
INTEGER          :: ii, jj, iVar, iGP, jGP
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

S = 0.0 ! S(1:nVar,nElemsX,nElemsY)

!int_{iixjj} S(1:nVar,ii,jj) dxdy

IF (GravitationalPotentialFlag .GT. 0) THEN

  SELECT CASE (ABS(Reconstruction))
    CASE(1,10,2,20,21,22,23,24,25,26,27,28,29)
      DO jj=1,nElemsY
        DO ii=1,nElemsX
          S_in_qp(1:nVar,nGPs,nGPs,ii,jj) = SourceFuncPrimitive( W(1:nVar,ii,jj) , MB(:,ii,jj) )
        END DO
      END DO
    CASE(3)
      PRINT*, "No option for reconstructing a specific kind of variable in source"
      STOP
      DO iVar=1,nVar
        DO jj=1,nElemsY
          DO ii=-nGhosts,nElemsX+nGhosts+1
             CALL WENO3_SecondSweep( W(iVar,ii,jj-nGhosts:jj+nGhosts) , Vtemp(iVar,ii,jj,1:nGPs) )
          END DO
        END DO
      END DO
      DO jj=1,nElemsY
        DO ii=1,nElemsX
          DO iGP=1,nGPs
            DO iVar=1,nVar
              CALL WENO3_SecondSweep( Vtemp(iVar,ii-nGhosts:ii+nGhosts,jj,iGP) , Vtemp2(iVar,1:nGPs,iGP,ii,jj) )
            END DO
            DO jGP=1,nGPs
              S_in_qp(1:nVar,jGP,iGP,ii,jj) = SourceFuncPrimitive( Vtemp2(1:nVar,jGP,iGP,ii,jj) , MGP(:,ii,jj,jGP,iGP)  )
            END DO 
          END DO 
        END DO
      END DO
    CASE(4,5)
      PRINT*, "No option for reconstructing a specific kind of variable in source"
      STOP
      DO iVar=1,nVar
        DO jj=1,nElemsY
          DO ii=-nGhosts,nElemsX+nGhosts+1
             CALL WENO5_SecondSweep( W(iVar,ii,jj-nGhosts:jj+nGhosts) , Vtemp(iVar,ii,jj,1:nGPs) )
          END DO
        END DO
      END DO
      DO jj=1,nElemsY
        DO ii=1,nElemsX
          DO iGP=1,nGPs
            DO iVar=1,nVar
              CALL WENO5_SecondSweep( Vtemp(iVar,ii-nGhosts:ii+nGhosts,jj,iGP) , Vtemp2(iVar,1:nGPs,iGP,ii,jj) )
            END DO
            DO jGP=1,nGPs
              S_in_qp(1:nVar,jGP,iGP,ii,jj) = SourceFuncPrimitive( Vtemp2(1:nVar,jGP,iGP,ii,jj) , MGP(:,ii,jj,jGP,iGP)  )
            END DO 
          END DO 
        END DO
      END DO
    CASE(7)
      PRINT*, "No option for reconstructing a specific kind of variable in source"
      STOP
      DO iVar=1,nVar
        DO jj=1,nElemsY
          DO ii=-nGhosts,nElemsX+nGhosts+1
             CALL WENO7_SecondSweep( W(iVar,ii,jj-nGhosts:jj+nGhosts) , Vtemp(iVar,ii,jj,1:nGPs) )
          END DO
        END DO
      END DO
      DO jj=1,nElemsY
        DO ii=1,nElemsX
          DO iGP=1,nGPs
            DO iVar=1,nVar
              CALL WENO7_SecondSweep( Vtemp(iVar,ii-nGhosts:ii+nGhosts,jj,iGP) , Vtemp2(iVar,1:nGPs,iGP,ii,jj) )
            END DO
            DO jGP=1,nGPs
              S_in_qp(1:nVar,jGP,iGP,ii,jj) = SourceFuncPrimitive( Vtemp2(1:nVar,jGP,iGP,ii,jj) , MGP(:,ii,jj,jGP,iGP)  )
            END DO 
          END DO 
        END DO
      END DO

    CASE DEFAULT
      ErrorMessage = "Reconstruction not implemented"
      WRITE(*,*) ErrorMessage
      STOP
  END SELECT


  DO jj=1,nElemsY
    DO ii=1,nElemsX
      DO iGP=1,nGPs
        DO jGP=1,nGPs
          S(1:nVar,ii,jj) = S(1:nVar,ii,jj) + WeightsGP(iGP,jGP) * S_in_qp(1:nVar,iGP,jGP,ii,jj) 
        END DO 
      END DO 
    END DO
  END DO

END IF
!-------------------------------------------------------------------------------!
END SUBROUTINE SourceTerm_Primitive_INPUT_PRIMITIVE
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
FUNCTION SourceFuncPrimitive(Q,X) RESULT(S) !*ALERT, it should be in conserved variables but I never tested it
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nVar
IMPLICIT NONE
REAL, DIMENSION(1:nVar), INTENT(IN)  :: Q 
REAL, DIMENSION(1:nDims) , INTENT(IN)  :: X
REAL, DIMENSION(1:nVar) :: S 

S(1) = 0.
#ifdef MOMENTUMINPRIMITIVEVARIABLES
S(2) = -Q(1)*Gravitational_Potential_X(X) 
S(3) = -Q(1)*Gravitational_Potential_Y(X) 
#else
S(2) = -Gravitational_Potential_X(X) 
S(3) = -Gravitational_Potential_Y(X) 
#endif
S(4) = 0.0

!-------------------------------------------------------------------------------!
END FUNCTION SourceFuncPrimitive
#endif
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE BoundaryConditions(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsType
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: V
#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
#endif
#ifdef ACTIVEFLUX
!*For staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: W_X

!*For staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
#endif
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
INTEGER            :: idx_vx, idx_vy
REAL               :: x0, xc, xt
REAL               :: Prim_in(1:nVar), Prim_out(1:nVar)
REAL               :: Cons_in(1:nVar), Cons_out(1:nVar)
REAL               :: ConsRefState1(1:nVar), ConsRefState2(1:nVar)
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

idx_vx = 2
idx_vy = 3

!------------------------------!
! Left Boundary Conditions     !
!------------------------------!
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nElemsX-nGhosts+ii,jj)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,-nGhosts+ii,jj) = WC(1:nVar,nElemsX-nGhosts+ii,jj)
#endif
      END DO
    END DO
  CASE(2) ! Transmissive
#if(1==1)
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        ! PRINT*, "CORRECT BCs"
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nGhosts-ii+1,jj)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,-nGhosts+ii,jj) = WC(1:nVar,nGhosts-ii+1,jj)
#endif
      END DO
    END DO
#elif(1==0)
    PRINT*, "DEBUG BCs, do not run like this"
    STOP
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,1,jj)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,-nGhosts+ii,jj) = WC(1:nVar,1,jj)
#endif
      END DO
    END DO
#elif(1==0)
    PRINT*, "DEBUG BCs, do not run like this"
    STOP
    Prim_in(1:nVar) = PrimRefState4(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = Cons_in(1:nVar)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,-nGhosts+ii,jj) = Prim_in(1:nVar)
#endif
      END DO
    END DO
#endif
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState4(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = Cons_in(1:nVar)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,-nGhosts+ii,jj) = Prim_in(1:nVar)
#endif
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState4(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = Cons_out(1:nVar)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,-nGhosts+ii,jj) = Prim_out(1:nVar)
#endif
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nGhosts-ii+1,jj)
        U(idx_vx,-nGhosts+ii,jj) =-U(idx_vx,nGhosts-ii+1,jj)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,-nGhosts+ii,jj) = WC(1:nVar,nGhosts-ii+1,jj)
        WC(idx_vx,-nGhosts+ii,jj) =-WC(idx_vx,nGhosts-ii+1,jj)
#endif
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition L not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

#ifdef ACTIVEFLUX
!*Staggering in X direction
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=0,nGhosts 
        W_X(1:nVar,-nGhosts+ii,jj) = W_X(1:nVar,nElemsX-nGhosts+ii,jj) !*NB: NO CHANGE IN INDICES
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        W_X(1:nVar,-nGhosts+ii,jj) = W_X(1:nVar,nGhosts-ii+1+1,jj) !*NB: +1 at the RHS
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState4(1:nVar)
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        W_X(1:nVar,-nGhosts+ii,jj) = Prim_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState4(1:nVar)
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        W_X(1:nVar,-nGhosts+ii,jj) = Prim_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        W_X(1:nVar,-nGhosts+ii,jj) = W_X(1:nVar,nGhosts-ii+1+1,jj) !*NB:+1 at the RHS
        W_X(idx_vx,-nGhosts+ii,jj) =-W_X(idx_vx,nGhosts-ii+1+1,jj) !*NB:+1 at the RHS
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition L not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!*Staggering in Y direction
!*NB: +1 in loop in Y direction
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=0,nGhosts
        W_Y(1:nVar,-nGhosts+ii,jj) = W_Y(1:nVar,nElemsX-nGhosts+ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=0,nGhosts
        W_Y(1:nVar,-nGhosts+ii,jj) = W_Y(1:nVar,nGhosts-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState4(1:nVar)
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=0,nGhosts
        W_Y(1:nVar,-nGhosts+ii,jj) = Prim_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState4(1:nVar)
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=0,nGhosts
        W_Y(1:nVar,-nGhosts+ii,jj) = Prim_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=0,nGhosts
        W_Y(1:nVar,-nGhosts+ii,jj) = W_Y(1:nVar,nGhosts-ii+1,jj)
        W_Y(idx_vx,-nGhosts+ii,jj) =-W_Y(idx_vx,nGhosts-ii+1,jj)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition L not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

#endif


!------------------------------!
! Right Boundary Conditions    !
!------------------------------!
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,ii,jj)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,nElemsX+ii,jj) = WC(1:nVar,ii,jj)
#endif
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY    
      DO ii=1,nGhosts+1
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,nElemsX-ii+1,jj)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,nElemsX+ii,jj) = WC(1:nVar,nElemsX-ii+1,jj)
#endif
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState2(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        U(1:nVar,nElemsX+ii,jj) = Cons_in(1:nVar)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,nElemsX+ii,jj) = Prim_in(1:nVar)
#endif
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState2(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        U(1:nVar,nElemsX+ii,jj) = Cons_out(1:nVar)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,nElemsX+ii,jj) = Prim_out(1:nVar)
#endif
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,nElemsX-ii+1,jj)
        U(idx_vx,nElemsX+ii,jj) =-U(idx_vx,nElemsX-ii+1,jj)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,nElemsX+ii,jj) = WC(1:nVar,nElemsX-ii+1,jj)
        WC(idx_vx,nElemsX+ii,jj) =-WC(idx_vx,nElemsX-ii+1,jj)
#endif
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition R not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

#ifdef ACTIVEFLUX

!*Staggering in X direction
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        W_X(1:nVar,nElemsX+ii+1,jj) = W_X(1:nVar,ii+1,jj) !*NB:+1 at both sides
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY    
      DO ii=1,nGhosts+1
        W_X(1:nVar,nElemsX+ii+1,jj) = W_X(1:nVar,nElemsX-ii+1,jj) !*NB:+1 at LHS
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState2(1:nVar)
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        W_X(1:nVar,nElemsX+ii+1,jj) = Prim_in(1:nVar) !*NB:+1 at LHS
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState2(1:nVar)
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        W_X(1:nVar,nElemsX+ii+1,jj) = Prim_out(1:nVar) !*NB:+1 at LHS
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        W_X(1:nVar,nElemsX+ii+1,jj) = W_X(1:nVar,nElemsX-ii+1,jj) !*NB:+1 at LHS
        W_X(idx_vx,nElemsX+ii+1,jj) =-W_X(idx_vx,nElemsX-ii+1,jj) !*NB:+1 at LHS
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition R not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!*Staggering in Y direction
!*NB: +1 in loop in Y direction
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=1,nGhosts+1
        W_Y(1:nVar,nElemsX+ii,jj) = W_Y(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction    
      DO ii=1,nGhosts+1
        W_Y(1:nVar,nElemsX+ii,jj) = W_Y(1:nVar,nElemsX-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState2(1:nVar)
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=1,nGhosts+1
        W_Y(1:nVar,nElemsX+ii,jj) = Prim_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState2(1:nVar)
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=1,nGhosts+1
        W_Y(1:nVar,nElemsX+ii,jj) = Prim_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=1,nGhosts+1
        W_Y(1:nVar,nElemsX+ii,jj) = W_Y(1:nVar,nElemsX-ii+1,jj)
        W_Y(idx_vx,nElemsX+ii,jj) =-W_Y(idx_vx,nElemsX-ii+1,jj)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition R not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

#endif

!------------------------------!
! Top Boundary Conditions      !
!------------------------------!
SELECT CASE(BoundaryConditionsType(3))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        U(1:nVar,ii,nElemsY+jj) = U(1:nVar,ii,jj)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,ii,nElemsY+jj) = WC(1:nVar,ii,jj)
#endif
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        U(1:nVar,ii,nElemsY+jj) = U(1:nVar,ii,nElemsY-jj+1)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,ii,nElemsY+jj) = WC(1:nVar,ii,nElemsY-jj+1)
#endif
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState3(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        U(1:nVar,ii,nElemsY+jj) = Cons_in(1:nVar)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,ii,nElemsY+jj) = Prim_in(1:nVar)
#endif
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState3(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        U(1:nVar,ii,nElemsY+jj) = Cons_out(1:nVar)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,ii,nElemsY+jj) = Prim_out(1:nVar)
#endif
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        U(1:nVar,ii,nElemsY+jj) = U(1:nVar,ii,nElemsY-jj+1)
        U(idx_vy,ii,nElemsY+jj) =-U(idx_vy,ii,nElemsY-jj+1)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,ii,nElemsY+jj) = WC(1:nVar,ii,nElemsY-jj+1)
        WC(idx_vy,ii,nElemsY+jj) =-WC(idx_vy,ii,nElemsY-jj+1)
#endif
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition T not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

#ifdef ACTIVEFLUX
!*Staggering in X direction
!*NB: +1 in loop in X direction
SELECT CASE(BoundaryConditionsType(3))
  CASE(1) ! Periodic
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=1,nGhosts+1
        W_X(1:nVar,ii,nElemsY+jj) = W_X(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=1,nGhosts+1
        W_X(1:nVar,ii,nElemsY+jj) = W_X(1:nVar,ii,nElemsY-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState3(1:nVar)
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=1,nGhosts+1
        W_X(1:nVar,ii,nElemsY+jj) = Prim_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState3(1:nVar)
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=1,nGhosts+1
        W_X(1:nVar,ii,nElemsY+jj) = Prim_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=1,nGhosts+1
        W_X(1:nVar,ii,nElemsY+jj) = W_X(1:nVar,ii,nElemsY-jj+1)
        W_X(idx_vy,ii,nElemsY+jj) =-W_X(idx_vy,ii,nElemsY-jj+1)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition T not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!*Staggering in Y direction
SELECT CASE(BoundaryConditionsType(3))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        W_Y(1:nVar,ii,nElemsY+jj+1) = W_Y(1:nVar,ii,jj+1) !*NB:+1 at both sides
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        W_Y(1:nVar,ii,nElemsY+jj+1) = W_Y(1:nVar,ii,nElemsY-jj+1) !*NB:+1 at LHS
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState3(1:nVar)
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        W_Y(1:nVar,ii,nElemsY+jj+1) = Prim_in(1:nVar) !*NB:+1 at LHS
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState3(1:nVar)
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        W_Y(1:nVar,ii,nElemsY+jj+1) = Prim_out(1:nVar) !*NB:+1 at LHS
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        W_Y(1:nVar,ii,nElemsY+jj+1) = W_Y(1:nVar,ii,nElemsY-jj+1) !"NB:+1 at LHS
        W_Y(idx_vy,ii,nElemsY+jj+1) =-W_Y(idx_vy,ii,nElemsY-jj+1) !"NB:+1 at LHS
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition T not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT
#endif


!------------------------------!
! Bottom Boundary Conditions   !
!------------------------------!
SELECT CASE(BoundaryConditionsType(1))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nElemsY-nGhosts+jj)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,ii,-nGhosts+jj) = WC(1:nVar,ii,nElemsY-nGhosts+jj)
#endif
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nGhosts-jj+1)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,ii,-nGhosts+jj) = WC(1:nVar,ii,nGhosts-jj+1)
#endif
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState1(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = Cons_in(1:nVar)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,ii,-nGhosts+jj) = Prim_in(1:nVar)
#endif
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState1(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = Cons_out(1:nVar)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,ii,-nGhosts+jj) = Prim_out(1:nVar)
#endif
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nGhosts-jj+1)
        U(idx_vy,ii,-nGhosts+jj) =-U(idx_vy,ii,nGhosts-jj+1)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,ii,-nGhosts+jj) = WC(1:nVar,ii,nGhosts-jj+1)
        WC(idx_vy,ii,-nGhosts+jj) =-WC(idx_vy,ii,nGhosts-jj+1)
#endif
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition B not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


#ifdef ACTIVEFLUX
!*Staggering in X direction
!*NB: +1 in loop in X direction
SELECT CASE(BoundaryConditionsType(1))
  CASE(1) ! Periodic
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=0,nGhosts
        W_X(1:nVar,ii,-nGhosts+jj) = W_X(1:nVar,ii,nElemsY-nGhosts+jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=0,nGhosts
        W_X(1:nVar,ii,-nGhosts+jj) = W_X(1:nVar,ii,nGhosts-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState1(1:nVar)
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=0,nGhosts
        W_X(1:nVar,ii,-nGhosts+jj) = Prim_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState1(1:nVar)
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=0,nGhosts
        W_X(1:nVar,ii,-nGhosts+jj) = Prim_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=0,nGhosts
        W_X(1:nVar,ii,-nGhosts+jj) = W_X(1:nVar,ii,nGhosts-jj+1)
        W_X(idx_vy,ii,-nGhosts+jj) =-W_X(idx_vy,ii,nGhosts-jj+1)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition B not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!*Staggering in Y direction
SELECT CASE(BoundaryConditionsType(1))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        W_Y(1:nVar,ii,-nGhosts+jj) = W_Y(1:nVar,ii,nElemsY-nGhosts+jj) !*NB: NO CHANGE IN INDICES
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        W_Y(1:nVar,ii,-nGhosts+jj) = W_Y(1:nVar,ii,nGhosts-jj+1+1) !*NB:+1 at the RHS
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState1(1:nVar)
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        W_Y(1:nVar,ii,-nGhosts+jj) = Prim_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState1(1:nVar)
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        W_Y(1:nVar,ii,-nGhosts+jj) = Prim_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        W_Y(1:nVar,ii,-nGhosts+jj) = W_Y(1:nVar,ii,nGhosts-jj+1+1) !*NB: +1 at the RHS
        W_Y(idx_vy,ii,-nGhosts+jj) =-W_Y(idx_vy,ii,nGhosts-jj+1+1) !*NB: +1 at the RHS
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition B not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

#endif

!------------------------------!
! Left Corners Boundary Conditions!
!------------------------------!
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    !*Bottom left
    DO jj=-nGhosts,0
      DO ii=0,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nElemsX-nGhosts+ii,jj)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,-nGhosts+ii,jj) = WC(1:nVar,nElemsX-nGhosts+ii,jj)
#endif
      END DO
    END DO

    !*Top left
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=0,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nElemsX-nGhosts+ii,jj)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,-nGhosts+ii,jj) = WC(1:nVar,nElemsX-nGhosts+ii,jj)
#endif
      END DO
    END DO
END SELECT

#ifdef ACTIVEFLUX
!*Staggering in X direction
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    !*Bottom left
    DO jj=-nGhosts,0
      DO ii=0,nGhosts
        W_X(1:nVar,-nGhosts+ii,jj) = W_X(1:nVar,nElemsX-nGhosts+ii,jj) !*NB: NO CHANGE IN INDICES
      END DO
    END DO

    !*Top left
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=0,nGhosts
        W_X(1:nVar,-nGhosts+ii,jj) = W_X(1:nVar,nElemsX-nGhosts+ii,jj) !*NB: NO CHANGE IN INDICES
      END DO
    END DO
END SELECT

!*Staggering in Y direction
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    !*Bottom left
    DO jj=-nGhosts,0
      DO ii=0,nGhosts
        W_Y(1:nVar,-nGhosts+ii,jj) = W_Y(1:nVar,nElemsX-nGhosts+ii,jj) !*NB: NO CHANGE IN INDICES
      END DO
    END DO

    !*Top left
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=0,nGhosts
        W_Y(1:nVar,-nGhosts+ii,jj+1) = W_Y(1:nVar,nElemsX-nGhosts+ii,jj+1) !*NB:+1 at both sides
      END DO
    END DO
END SELECT


#endif


!------------------------------!
! Right Corners Boundary Conditions    !
!------------------------------!
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    !*Bottom right
    DO jj=-nGhosts,0
      DO ii=1,nGhosts+1
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,ii,jj)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,nElemsX+ii,jj) = WC(1:nVar,ii,jj)
#endif
      END DO
    END DO

    !*Top right
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=1,nGhosts+1
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,ii,jj)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,nElemsX+ii,jj) = WC(1:nVar,ii,jj)
#endif
      END DO
    END DO
END SELECT

#ifdef ACTIVEFLUX
!*Staggering in X direction
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    !*Bottom right
    DO jj=-nGhosts,0
      DO ii=1,nGhosts+1
        W_X(1:nVar,nElemsX+ii+1,jj) = W_X(1:nVar,ii+1,jj) !*NB:+1 at both sides
      END DO
    END DO

    !*Top right
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=1,nGhosts+1
        W_X(1:nVar,nElemsX+ii+1,jj) = W_X(1:nVar,ii+1,jj) !*NB:+1 at both sides
      END DO
    END DO
END SELECT

!*Staggering in Y direction
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    !*Bottom right
    DO jj=-nGhosts,0
      DO ii=1,nGhosts+1
        W_Y(1:nVar,nElemsX+ii,jj) = W_Y(1:nVar,ii,jj) !*NB: NO CHANGE IN INDICES
      END DO
    END DO

    !*Top right
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=1,nGhosts+1
        W_Y(1:nVar,nElemsX+ii,jj+1) = W_Y(1:nVar,ii,jj+1) !*NB:+1 at both sides
      END DO
    END DO
END SELECT


#endif

#if(1==0)
PRINT*, "It did not solve the problem with left boundary of tests -302 and -304."
PRINT*, "Stop in equation in BoundaryConditions"
STOP

IF (InitialCondition .EQ. -304) THEN
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,1,jj)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,-nGhosts+ii,jj) = WC(1:nVar,1,jj)
#endif
      END DO
    END DO
    DO jj=1,nElemsY    
      DO ii=1,nGhosts+1
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,nElemsX,jj)
#ifdef CENTEREDPRIMITIVE
        WC(1:nVar,nElemsX+ii,jj) = WC(1:nVar,nElemsX,jj)
#endif
      END DO
    END DO
END IF
#endif


DO jj=1,nElemsY
  DO ii=1,nElemsX
    CALL ConsToPrim(U(1:nVar,ii,jj),V(1:nVar,ii,jj))
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE BoundaryConditions
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE BoundaryConditions_On_Cell_Centered_Data_Primitive_INPUT(W_in_out,t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsType
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: t
REAL,INTENT(INOUT) :: W_in_out(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) !*Like standard, cell centered
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
INTEGER            :: idx_vx, idx_vy
REAL               :: x0, xc, xt
REAL               :: Prim_in(1:nVar), Prim_out(1:nVar)
REAL               :: Cons_in(1:nVar), Cons_out(1:nVar)
REAL               :: ConsRefState1(1:nVar), ConsRefState2(1:nVar)
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

idx_vx = 2
idx_vy = 3

!------------------------------!
! Left Boundary Conditions     !
!------------------------------!
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        W_in_out(1:nVar,-nGhosts+ii,jj) = W_in_out(1:nVar,nElemsX-nGhosts+ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        W_in_out(1:nVar,-nGhosts+ii,jj) = W_in_out(1:nVar,nGhosts-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState4(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        W_in_out(1:nVar,-nGhosts+ii,jj) = Prim_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState4(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        W_in_out(1:nVar,-nGhosts+ii,jj) = Prim_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        W_in_out(1:nVar,-nGhosts+ii,jj) = W_in_out(1:nVar,nGhosts-ii+1,jj)
        W_in_out(idx_vx,-nGhosts+ii,jj) =-W_in_out(idx_vx,nGhosts-ii+1,jj)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition L not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!------------------------------!
! Right Boundary Conditions    !
!------------------------------!
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        W_in_out(1:nVar,nElemsX+ii,jj) = W_in_out(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY    
      DO ii=1,nGhosts+1
        W_in_out(1:nVar,nElemsX+ii,jj) = W_in_out(1:nVar,nElemsX-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState2(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        W_in_out(1:nVar,nElemsX+ii,jj) = Prim_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState2(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        W_in_out(1:nVar,nElemsX+ii,jj) = Prim_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        W_in_out(1:nVar,nElemsX+ii,jj) = W_in_out(1:nVar,nElemsX-ii+1,jj)
        W_in_out(idx_vx,nElemsX+ii,jj) =-W_in_out(idx_vx,nElemsX-ii+1,jj)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition R not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT



!------------------------------!
! Top Boundary Conditions      !
!------------------------------!
SELECT CASE(BoundaryConditionsType(3))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        W_in_out(1:nVar,ii,nElemsY+jj) = W_in_out(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        W_in_out(1:nVar,ii,nElemsY+jj) = W_in_out(1:nVar,ii,nElemsY-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState3(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        W_in_out(1:nVar,ii,nElemsY+jj) = Prim_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState3(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        W_in_out(1:nVar,ii,nElemsY+jj) = Prim_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        W_in_out(1:nVar,ii,nElemsY+jj) = W_in_out(1:nVar,ii,nElemsY-jj+1)
        W_in_out(idx_vy,ii,nElemsY+jj) =-W_in_out(idx_vy,ii,nElemsY-jj+1)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition T not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT



!------------------------------!
! Bottom Boundary Conditions   !
!------------------------------!
SELECT CASE(BoundaryConditionsType(1))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        W_in_out(1:nVar,ii,-nGhosts+jj) = W_in_out(1:nVar,ii,nElemsY-nGhosts+jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        W_in_out(1:nVar,ii,-nGhosts+jj) = W_in_out(1:nVar,ii,nGhosts-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState1(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        W_in_out(1:nVar,ii,-nGhosts+jj) = Prim_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState1(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        W_in_out(1:nVar,ii,-nGhosts+jj) = Prim_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        W_in_out(1:nVar,ii,-nGhosts+jj) = W_in_out(1:nVar,ii,nGhosts-jj+1)
        W_in_out(idx_vy,ii,-nGhosts+jj) =-W_in_out(idx_vy,ii,nGhosts-jj+1)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition B not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT



!------------------------------!
! Left Corners Boundary Conditions!
!------------------------------!
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    !*Bottom left
    DO jj=-nGhosts,0
      DO ii=0,nGhosts
        W_in_out(1:nVar,-nGhosts+ii,jj) = W_in_out(1:nVar,nElemsX-nGhosts+ii,jj)
      END DO
    END DO

    !*Top left
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=0,nGhosts
        W_in_out(1:nVar,-nGhosts+ii,jj) = W_in_out(1:nVar,nElemsX-nGhosts+ii,jj)
      END DO
    END DO
END SELECT



!------------------------------!
! Right Corners Boundary Conditions    !
!------------------------------!
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    !*Bottom right
    DO jj=-nGhosts,0
      DO ii=1,nGhosts+1
        W_in_out(1:nVar,nElemsX+ii,jj) = W_in_out(1:nVar,ii,jj)
      END DO
    END DO

    !*Top right
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=1,nGhosts+1
        W_in_out(1:nVar,nElemsX+ii,jj) = W_in_out(1:nVar,ii,jj)
      END DO
    END DO
END SELECT


!-------------------------------------------------------------------------------!
END SUBROUTINE BoundaryConditions_On_Cell_Centered_Data_Primitive_INPUT
!===============================================================================!
!
!
!
!===============================================================================!
FUNCTION TimeStep(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: CFL
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: LambdaMaxX
USE MOD_FiniteVolume2D_vars,ONLY: LambdaMaxY
USE MOD_FiniteVolume2D_vars,ONLY: MIN_TIMESTEP
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: NormVectX
USE MOD_FiniteVolume2D_vars,ONLY: TangVectX
USE MOD_FiniteVolume2D_vars,ONLY: NormVectY
USE MOD_FiniteVolume2D_vars,ONLY: TangVectY
USE MOD_FiniteVolume2D_vars,ONLY: WhichSpeedEstimateForDtComputation
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_Reconstruction,      ONLY: ReconstructionX_INPUT_CONSERVED
USE MOD_Reconstruction,      ONLY: ReconstructionY_INPUT_CONSERVED
#ifdef ACTIVEFLUX
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
#ifdef IMEX
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: min_p
#endif
#endif
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
USE exact_riemann_mod      ,ONLY: exact_riemann
#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, INTENT(IN)    :: t
REAL                :: TimeStep
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL    :: FastestWaveXInAbsoluteValue, FastestWaveYInAbsoluteValue
REAL    :: c, ro, vx, vy, p
REAL    :: Prim(1:nVar)
INTEGER :: ii, jj
REAL    :: PrimL(1:nVar), PrimR(1:nVar)
REAL    :: PrimLL(1:nVar), PrimRR(1:nVar)
REAL    :: NormVect(nDims), TangVect(nDims)
#ifdef BLENDINGDT
REAL    :: dt_explicit, theta
#endif
!-------------------------------------------------------------------------------!

LambdaMaxX = 0.0
LambdaMaxY = 0.0
TimeStep = HUGE(1.0)



#ifdef ACTIVEFLUX
!*Staggering X
DO jj=1,nElemsY
  DO ii=1,nElemsX+1 !*NB:+1 in X direction
    Prim(1:nVar)=W_X(1:nVar,ii,jj)

    !*NB: This should be ok also for defined(RELAXATION) && !(defined(IMEX)) 
    CALL WaveSpeeds2D(Prim(1:nVar),fastest_in_absolute_valuex=FastestWaveXInAbsoluteValue,fastest_in_absolute_valuey=FastestWaveYInAbsoluteValue)

    LambdaMaxX = MAX(LambdaMaxX,ABS(FastestWaveXInAbsoluteValue))
    LambdaMaxY = MAX(LambdaMaxY,ABS(FastestWaveYInAbsoluteValue))

    TimeStep  = MIN(TimeStep,MESH_DX(1)/LambdaMaxX,MESH_DX(2)/LambdaMaxY)
  END DO
END DO


!*Staggering Y
DO jj=1,nElemsY+1 !*NB:+1 in Y direction
  DO ii=1,nElemsX
    Prim(1:nVar)=W_Y(1:nVar,ii,jj)

    !*NB: This should be ok also for defined(RELAXATION) && !(defined(IMEX)) 
    CALL WaveSpeeds2D(Prim(1:nVar),fastest_in_absolute_valuex=FastestWaveXInAbsoluteValue,fastest_in_absolute_valuey=FastestWaveYInAbsoluteValue)

    LambdaMaxX = MAX(LambdaMaxX,ABS(FastestWaveXInAbsoluteValue))
    LambdaMaxY = MAX(LambdaMaxY,ABS(FastestWaveYInAbsoluteValue))

    TimeStep  = MIN(TimeStep,MESH_DX(1)/LambdaMaxX,MESH_DX(2)/LambdaMaxY)
  END DO
END DO

#elif defined(CENTEREDPRIMITIVE)

DO jj=1,nElemsY
  DO ii=1,nElemsX
    Prim(1:nVar)=WC(1:nVar,ii,jj)

    !*NB: This should be ok also for defined(RELAXATION) && !(defined(IMEX)) 
    CALL WaveSpeeds2D(Prim(1:nVar),fastest_in_absolute_valuex=FastestWaveXInAbsoluteValue,fastest_in_absolute_valuey=FastestWaveYInAbsoluteValue)

    LambdaMaxX = MAX(LambdaMaxX,ABS(FastestWaveXInAbsoluteValue))
    LambdaMaxY = MAX(LambdaMaxY,ABS(FastestWaveYInAbsoluteValue))

    TimeStep  = MIN(TimeStep,MESH_DX(1)/LambdaMaxX,MESH_DX(2)/LambdaMaxY)
  END DO
END DO

#else

!*=============================
!*IF NOT ACTIVEFLUX, THEN CONSIDER THE CONSERVED VARIABLES
!*=============================

#ifdef MULTIFLUID
PRINT*, "You should not be here in TimeStep"
STOP
#endif

#if(1==1)
DO jj=1,nElemsY
  DO ii=1,nElemsX
    CALL ConsToPrim(U(1:nVar,ii,jj),Prim(1:nVar))

    !*NB: This should be ok also for defined(RELAXATION) && !(defined(IMEX)) 
    CALL WaveSpeeds2D(Prim(1:nVar),fastest_in_absolute_valuex=FastestWaveXInAbsoluteValue,fastest_in_absolute_valuey=FastestWaveYInAbsoluteValue)

    LambdaMaxX = MAX(LambdaMaxX,ABS(FastestWaveXInAbsoluteValue))
    LambdaMaxY = MAX(LambdaMaxY,ABS(FastestWaveYInAbsoluteValue))
    TimeStep  = MIN(TimeStep,MESH_DX(1)/LambdaMaxX,MESH_DX(2)/LambdaMaxY)
  END DO
END DO
#endif

#if(1==0)
PRINT*, "I DO NOT WANT TO RUN IN THIS MODE (TIME STEP COMPUTATION AS ALEX)"
PRINT*, "Stopping in timestep in equation"
STOP
CALL BoundaryConditions(t)

CALL ReconstructionX_INPUT_CONSERVED(U( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&    !*Y bounds vector  
                          & 0,nElemsX+1, &                 !*X bounds loop  
                          & 1,nElemsY)                   !*Y bounds loop    

DO jj=1,nElemsY
  DO ii=0,nElemsX
    CALL ConsToPrim(WP(1:nVar,1,ii+0,jj),Prim(1:nVar))
    !*NB: This should be ok also for defined(RELAXATION) && !(defined(IMEX)) 
    CALL WaveSpeeds2D(Prim(1:nVar),fastest_in_absolute_valuex=FastestWaveXInAbsoluteValue,fastest_in_absolute_valuey=FastestWaveYInAbsoluteValue)

    LambdaMaxX = MAX(LambdaMaxX,ABS(FastestWaveXInAbsoluteValue))
    CALL ConsToPrim(WM(1:nVar,1,ii+1,jj),Prim(1:nVar))
    !*NB: This should be ok also for defined(RELAXATION) && !(defined(IMEX)) 
    CALL WaveSpeeds2D(Prim(1:nVar),fastest_in_absolute_valuex=FastestWaveXInAbsoluteValue,fastest_in_absolute_valuey=FastestWaveYInAbsoluteValue)

    LambdaMaxX = MAX(LambdaMaxX,ABS(FastestWaveXInAbsoluteValue))
  END DO
END DO



CALL ReconstructionY_INPUT_CONSERVED(U( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&  !*Y bounds vector 
                          & 1,nElemsX, &                   !*X bounds loop
                          & 0,nElemsY+1)                 !*Y bounds loop    

DO jj=0,nElemsY
  DO ii=1,nElemsX
    CALL ConsToPrim(WP(1:nVar,1,ii,jj+0),Prim(1:nVar))
    !*NB: This should be ok also for defined(RELAXATION) && !(defined(IMEX)) 
    CALL WaveSpeeds2D(Prim(1:nVar),fastest_in_absolute_valuex=FastestWaveXInAbsoluteValue,fastest_in_absolute_valuey=FastestWaveYInAbsoluteValue)

    LambdaMaxY = MAX(LambdaMaxY,ABS(FastestWaveYInAbsoluteValue))
    CALL ConsToPrim(WM(1:nVar,1,ii,jj+1),Prim(1:nVar))
    !*NB: This should be ok also for defined(RELAXATION) && !(defined(IMEX)) 
    CALL WaveSpeeds2D(Prim(1:nVar),fastest_in_absolute_valuex=FastestWaveXInAbsoluteValue,fastest_in_absolute_valuey=FastestWaveYInAbsoluteValue)

    LambdaMaxY = MAX(LambdaMaxY,ABS(FastestWaveYInAbsoluteValue))
  END DO
END DO

TimeStep  = MIN(TimeStep,MESH_DX(1)/LambdaMaxX,MESH_DX(2)/LambdaMaxY)
#endif


#endif


IF (WhichSpeedEstimateForDtComputation .EQ. 1) THEN


CALL BoundaryConditions(t)


#ifdef ACTIVEFLUX

!*W_X X
DO jj=1,nElemsY
  DO ii=0,nElemsX+1
    PrimL(1:nVar)=W_X(1:nVar,ii+0,jj)
    PrimR(1:nVar)=W_X(1:nVar,ii+1,jj)

    NormVect(1:nDims)=NormVectX(1:nDims)
    TangVect(1:nDims)=TangVectX(1:nDims)


    !*Rotate, actually here there is no rotation
    PrimLL(1) = PrimL(1)
    PrimLL(2) = NormVect(1)*PrimL(2) + NormVect(2)*PrimL(3)
    PrimLL(3) = TangVect(1)*PrimL(2) + TangVect(2)*PrimL(3)
    PrimLL(4) = PrimL(4)

    PrimRR(1) = PrimR(1)
    PrimRR(2) = NormVect(1)*PrimR(2) + NormVect(2)*PrimR(3)
    PrimRR(3) = TangVect(1)*PrimR(2) + TangVect(2)*PrimR(3)
    PrimRR(4) = PrimR(4)


    CALL exact_riemann(Gmm, &
                    & PrimLL(1),PrimLL(2),PrimLL(4),&
                    & PrimRR(1),PrimRR(2),PrimRR(4),&
                    & Prim(1),Prim(2),Prim(4))

    IF (Prim(2) .GT. 0.0) THEN
        Prim(3)=PrimLL(3)
    ELSE
        Prim(3)=PrimRR(3)
    ENDIF

    !*NB: This should be ok also for defined(RELAXATION) && !(defined(IMEX)) 
    CALL WaveSpeeds2D(Prim(1:nVar),fastest_in_absolute_valuex=FastestWaveXInAbsoluteValue,fastest_in_absolute_valuey=FastestWaveYInAbsoluteValue)

    LambdaMaxX = MAX(LambdaMaxX,ABS(FastestWaveXInAbsoluteValue))
    LambdaMaxY = MAX(LambdaMaxY,ABS(FastestWaveYInAbsoluteValue))

    TimeStep  = MIN(TimeStep,MESH_DX(1)/LambdaMaxX,MESH_DX(2)/LambdaMaxY)

  END DO
END DO

!*W_X Y
DO jj=0,nElemsY
  DO ii=1,nElemsX+1
    PrimL(1:nVar)=W_X(1:nVar,ii,jj+0)
    PrimR(1:nVar)=W_X(1:nVar,ii,jj+1)

    NormVect(1:nDims)=NormVectY(1:nDims)
    TangVect(1:nDims)=TangVectY(1:nDims)

    !*Rotate, actually here there is no rotation
    PrimLL(1) = PrimL(1)
    PrimLL(2) = NormVect(1)*PrimL(2) + NormVect(2)*PrimL(3)
    PrimLL(3) = TangVect(1)*PrimL(2) + TangVect(2)*PrimL(3)
    PrimLL(4) = PrimL(4)

    PrimRR(1) = PrimR(1)
    PrimRR(2) = NormVect(1)*PrimR(2) + NormVect(2)*PrimR(3)
    PrimRR(3) = TangVect(1)*PrimR(2) + TangVect(2)*PrimR(3)
    PrimRR(4) = PrimR(4)

    CALL exact_riemann(Gmm, &
                    & PrimLL(1),PrimLL(2),PrimLL(4),&
                    & PrimRR(1),PrimRR(2),PrimRR(4),&
                    & Prim(1),Prim(2),Prim(4))

    IF (Prim(2) .GT. 0.0) THEN
        Prim(3)=PrimLL(3)
    ELSE
        Prim(3)=PrimRR(3)
    ENDIF

    !*NB: This should be ok also for defined(RELAXATION) && !(defined(IMEX)) 
    CALL WaveSpeeds2D(Prim(1:nVar),fastest_in_absolute_valuex=FastestWaveXInAbsoluteValue,fastest_in_absolute_valuey=FastestWaveYInAbsoluteValue)

    LambdaMaxX = MAX(LambdaMaxX,ABS(FastestWaveXInAbsoluteValue))
    LambdaMaxY = MAX(LambdaMaxY,ABS(FastestWaveYInAbsoluteValue))

    TimeStep  = MIN(TimeStep,MESH_DX(1)/LambdaMaxX,MESH_DX(2)/LambdaMaxY)

  END DO
END DO


!*====

!*W_Y X
DO jj=1,nElemsY+1
  DO ii=0,nElemsX
    PrimL(1:nVar)=W_Y(1:nVar,ii+0,jj)
    PrimR(1:nVar)=W_Y(1:nVar,ii+1,jj)

    NormVect(1:nDims)=NormVectX(1:nDims)
    TangVect(1:nDims)=TangVectX(1:nDims)


    !*Rotate, actually here there is no rotation
    PrimLL(1) = PrimL(1)
    PrimLL(2) = NormVect(1)*PrimL(2) + NormVect(2)*PrimL(3)
    PrimLL(3) = TangVect(1)*PrimL(2) + TangVect(2)*PrimL(3)
    PrimLL(4) = PrimL(4)

    PrimRR(1) = PrimR(1)
    PrimRR(2) = NormVect(1)*PrimR(2) + NormVect(2)*PrimR(3)
    PrimRR(3) = TangVect(1)*PrimR(2) + TangVect(2)*PrimR(3)
    PrimRR(4) = PrimR(4)


    CALL exact_riemann(Gmm, &
                    & PrimLL(1),PrimLL(2),PrimLL(4),&
                    & PrimRR(1),PrimRR(2),PrimRR(4),&
                    & Prim(1),Prim(2),Prim(4))

    IF (Prim(2) .GT. 0.0) THEN
        Prim(3)=PrimLL(3)
    ELSE
        Prim(3)=PrimRR(3)
    ENDIF

    !*NB: This should be ok also for defined(RELAXATION) && !(defined(IMEX)) 
    CALL WaveSpeeds2D(Prim(1:nVar),fastest_in_absolute_valuex=FastestWaveXInAbsoluteValue,fastest_in_absolute_valuey=FastestWaveYInAbsoluteValue)

    LambdaMaxX = MAX(LambdaMaxX,ABS(FastestWaveXInAbsoluteValue))
    LambdaMaxY = MAX(LambdaMaxY,ABS(FastestWaveYInAbsoluteValue))

    TimeStep  = MIN(TimeStep,MESH_DX(1)/LambdaMaxX,MESH_DX(2)/LambdaMaxY)

  END DO
END DO

!*W_X Y
DO jj=0,nElemsY+1
  DO ii=1,nElemsX
    PrimL(1:nVar)=W_Y(1:nVar,ii,jj+0)
    PrimR(1:nVar)=W_Y(1:nVar,ii,jj+1)

    NormVect(1:nDims)=NormVectY(1:nDims)
    TangVect(1:nDims)=TangVectY(1:nDims)

    !*Rotate, actually here there is no rotation
    PrimLL(1) = PrimL(1)
    PrimLL(2) = NormVect(1)*PrimL(2) + NormVect(2)*PrimL(3)
    PrimLL(3) = TangVect(1)*PrimL(2) + TangVect(2)*PrimL(3)
    PrimLL(4) = PrimL(4)

    PrimRR(1) = PrimR(1)
    PrimRR(2) = NormVect(1)*PrimR(2) + NormVect(2)*PrimR(3)
    PrimRR(3) = TangVect(1)*PrimR(2) + TangVect(2)*PrimR(3)
    PrimRR(4) = PrimR(4)

    CALL exact_riemann(Gmm, &
                    & PrimLL(1),PrimLL(2),PrimLL(4),&
                    & PrimRR(1),PrimRR(2),PrimRR(4),&
                    & Prim(1),Prim(2),Prim(4))

    IF (Prim(2) .GT. 0.0) THEN
        Prim(3)=PrimLL(3)
    ELSE
        Prim(3)=PrimRR(3)
    ENDIF

    !*NB: This should be ok also for defined(RELAXATION) && !(defined(IMEX)) 
    CALL WaveSpeeds2D(Prim(1:nVar),fastest_in_absolute_valuex=FastestWaveXInAbsoluteValue,fastest_in_absolute_valuey=FastestWaveYInAbsoluteValue)

    LambdaMaxX = MAX(LambdaMaxX,ABS(FastestWaveXInAbsoluteValue))
    LambdaMaxY = MAX(LambdaMaxY,ABS(FastestWaveYInAbsoluteValue))

    TimeStep  = MIN(TimeStep,MESH_DX(1)/LambdaMaxX,MESH_DX(2)/LambdaMaxY)

  END DO
END DO

#elif defined(CENTEREDPRIMITIVE)


!*WC X
DO jj=1,nElemsY
  DO ii=0,nElemsX
    PrimL(1:nVar)=WC(1:nVar,ii+0,jj)
    PrimR(1:nVar)=WC(1:nVar,ii+1,jj)

    NormVect(1:nDims)=NormVectX(1:nDims)
    TangVect(1:nDims)=TangVectX(1:nDims)


    !*Rotate, actually here there is no rotation
    PrimLL(1) = PrimL(1)
    PrimLL(2) = NormVect(1)*PrimL(2) + NormVect(2)*PrimL(3)
    PrimLL(3) = TangVect(1)*PrimL(2) + TangVect(2)*PrimL(3)
    PrimLL(4) = PrimL(4)

    PrimRR(1) = PrimR(1)
    PrimRR(2) = NormVect(1)*PrimR(2) + NormVect(2)*PrimR(3)
    PrimRR(3) = TangVect(1)*PrimR(2) + TangVect(2)*PrimR(3)
    PrimRR(4) = PrimR(4)


    CALL exact_riemann(Gmm, &
                    & PrimLL(1),PrimLL(2),PrimLL(4),&
                    & PrimRR(1),PrimRR(2),PrimRR(4),&
                    & Prim(1),Prim(2),Prim(4))

    IF (Prim(2) .GT. 0.0) THEN
        Prim(3)=PrimLL(3)
    ELSE
        Prim(3)=PrimRR(3)
    ENDIF

    !*NB: This should be ok also for defined(RELAXATION) && !(defined(IMEX)) 
    CALL WaveSpeeds2D(Prim(1:nVar),fastest_in_absolute_valuex=FastestWaveXInAbsoluteValue,fastest_in_absolute_valuey=FastestWaveYInAbsoluteValue)

    LambdaMaxX = MAX(LambdaMaxX,ABS(FastestWaveXInAbsoluteValue))
    LambdaMaxY = MAX(LambdaMaxY,ABS(FastestWaveYInAbsoluteValue))

    TimeStep  = MIN(TimeStep,MESH_DX(1)/LambdaMaxX,MESH_DX(2)/LambdaMaxY)

  END DO
END DO

!*WC Y
DO jj=0,nElemsY
  DO ii=1,nElemsX
    PrimL(1:nVar)=WC(1:nVar,ii,jj+0)
    PrimR(1:nVar)=WC(1:nVar,ii,jj+1)

    NormVect(1:nDims)=NormVectY(1:nDims)
    TangVect(1:nDims)=TangVectY(1:nDims)

    !*Rotate, actually here there is no rotation
    PrimLL(1) = PrimL(1)
    PrimLL(2) = NormVect(1)*PrimL(2) + NormVect(2)*PrimL(3)
    PrimLL(3) = TangVect(1)*PrimL(2) + TangVect(2)*PrimL(3)
    PrimLL(4) = PrimL(4)

    PrimRR(1) = PrimR(1)
    PrimRR(2) = NormVect(1)*PrimR(2) + NormVect(2)*PrimR(3)
    PrimRR(3) = TangVect(1)*PrimR(2) + TangVect(2)*PrimR(3)
    PrimRR(4) = PrimR(4)

    CALL exact_riemann(Gmm, &
                    & PrimLL(1),PrimLL(2),PrimLL(4),&
                    & PrimRR(1),PrimRR(2),PrimRR(4),&
                    & Prim(1),Prim(2),Prim(4))

    IF (Prim(2) .GT. 0.0) THEN
        Prim(3)=PrimLL(3)
    ELSE
        Prim(3)=PrimRR(3)
    ENDIF

    !*NB: This should be ok also for defined(RELAXATION) && !(defined(IMEX)) 
    CALL WaveSpeeds2D(Prim(1:nVar),fastest_in_absolute_valuex=FastestWaveXInAbsoluteValue,fastest_in_absolute_valuey=FastestWaveYInAbsoluteValue)

    LambdaMaxX = MAX(LambdaMaxX,ABS(FastestWaveXInAbsoluteValue))
    LambdaMaxY = MAX(LambdaMaxY,ABS(FastestWaveYInAbsoluteValue))

    TimeStep  = MIN(TimeStep,MESH_DX(1)/LambdaMaxX,MESH_DX(2)/LambdaMaxY)

  END DO
END DO




#else


#ifdef MULTIFLUID
PRINT*, "You should not be here in TimeStep"
STOP
#endif

!*X
DO jj=1,nElemsY
  DO ii=0,nElemsX
    CALL ConsToPrim(U(1:nVar,ii+0,jj),PrimL(1:nVar))
    CALL ConsToPrim(U(1:nVar,ii+1,jj),PrimR(1:nVar))

    NormVect(1:nDims)=NormVectX(1:nDims)
    TangVect(1:nDims)=TangVectX(1:nDims)


    !*Rotate, actually here there is no rotation
    PrimLL(1) = PrimL(1)
    PrimLL(2) = NormVect(1)*PrimL(2) + NormVect(2)*PrimL(3)
    PrimLL(3) = TangVect(1)*PrimL(2) + TangVect(2)*PrimL(3)
    PrimLL(4) = PrimL(4)

    PrimRR(1) = PrimR(1)
    PrimRR(2) = NormVect(1)*PrimR(2) + NormVect(2)*PrimR(3)
    PrimRR(3) = TangVect(1)*PrimR(2) + TangVect(2)*PrimR(3)
    PrimRR(4) = PrimR(4)


    CALL exact_riemann(Gmm, &
                    & PrimLL(1),PrimLL(2),PrimLL(4),&
                    & PrimRR(1),PrimRR(2),PrimRR(4),&
                    & Prim(1),Prim(2),Prim(4))

    IF (Prim(2) .GT. 0.0) THEN
        Prim(3)=PrimLL(3)
    ELSE
        Prim(3)=PrimRR(3)
    ENDIF

    !*NB: This should be ok also for defined(RELAXATION) && !(defined(IMEX)) 
    CALL WaveSpeeds2D(Prim(1:nVar),fastest_in_absolute_valuex=FastestWaveXInAbsoluteValue,fastest_in_absolute_valuey=FastestWaveYInAbsoluteValue)

    LambdaMaxX = MAX(LambdaMaxX,ABS(FastestWaveXInAbsoluteValue))
    LambdaMaxY = MAX(LambdaMaxY,ABS(FastestWaveYInAbsoluteValue))

    TimeStep  = MIN(TimeStep,MESH_DX(1)/LambdaMaxX,MESH_DX(2)/LambdaMaxY)

  END DO
END DO

!*Y
DO jj=0,nElemsY
  DO ii=1,nElemsX
    CALL ConsToPrim(U(1:nVar,ii,jj+0),PrimL(1:nVar))
    CALL ConsToPrim(U(1:nVar,ii,jj+1),PrimR(1:nVar))

    NormVect(1:nDims)=NormVectY(1:nDims)
    TangVect(1:nDims)=TangVectY(1:nDims)

    !*Rotate, actually here there is no rotation
    PrimLL(1) = PrimL(1)
    PrimLL(2) = NormVect(1)*PrimL(2) + NormVect(2)*PrimL(3)
    PrimLL(3) = TangVect(1)*PrimL(2) + TangVect(2)*PrimL(3)
    PrimLL(4) = PrimL(4)

    PrimRR(1) = PrimR(1)
    PrimRR(2) = NormVect(1)*PrimR(2) + NormVect(2)*PrimR(3)
    PrimRR(3) = TangVect(1)*PrimR(2) + TangVect(2)*PrimR(3)
    PrimRR(4) = PrimR(4)

    CALL exact_riemann(Gmm, &
                    & PrimLL(1),PrimLL(2),PrimLL(4),&
                    & PrimRR(1),PrimRR(2),PrimRR(4),&
                    & Prim(1),Prim(2),Prim(4))

    IF (Prim(2) .GT. 0.0) THEN
        Prim(3)=PrimLL(3)
    ELSE
        Prim(3)=PrimRR(3)
    ENDIF

    !*NB: This should be ok also for defined(RELAXATION) && !(defined(IMEX)) 
    CALL WaveSpeeds2D(Prim(1:nVar),fastest_in_absolute_valuex=FastestWaveXInAbsoluteValue,fastest_in_absolute_valuey=FastestWaveYInAbsoluteValue)

    LambdaMaxX = MAX(LambdaMaxX,ABS(FastestWaveXInAbsoluteValue))
    LambdaMaxY = MAX(LambdaMaxY,ABS(FastestWaveYInAbsoluteValue))

    TimeStep  = MIN(TimeStep,MESH_DX(1)/LambdaMaxX,MESH_DX(2)/LambdaMaxY)

  END DO
END DO

#endif

END IF

#ifdef IMEX
#ifdef BLENDINGDT
dt_explicit=HUGE(1.0)
DO jj=1,nElemsY
  DO ii=1,nElemsX
    Prim(1:nVar)=WC(1:nVar,ii,jj)

    !*NB: This should be ok also for defined(RELAXATION) && !(defined(IMEX)) 
    CALL WaveSpeeds2D(Prim(1:nVar),fastestx_without_splitting=FastestWaveXInAbsoluteValue,fastesty_without_splitting=FastestWaveYInAbsoluteValue)

    LambdaMaxX = MAX(LambdaMaxX,ABS(FastestWaveXInAbsoluteValue))
    LambdaMaxY = MAX(LambdaMaxY,ABS(FastestWaveYInAbsoluteValue))
    dt_explicit  = MIN(dt_explicit,MESH_DX(1)/LambdaMaxX,MESH_DX(2)/LambdaMaxY)
  END DO
END DO

#ifdef RELAXATION
theta=Switching_Function(EPS_LM)
IF (TimeStep .GT. MESH_DX(1)*1000.0) THEN !*IF THE IMPLICIT DOES NOT MAKE SENSE
  TimeStep=dt_explicit
ELSE 
  TimeStep=(1.0-theta)*dt_explicit+theta*TimeStep
END IF
#else
TimeStep=dt_explicit
#endif

#endif
#endif

TimeStep = CFL*TimeStep

IF (TimeStep .LT. MIN_TIMESTEP) THEN
  TimeStep = MIN_TIMESTEP
END IF

!-------------------------------------------------------------------------------!
END FUNCTION TimeStep
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WaveSpeeds1D(Prim,slowest,fastest,fastest_in_absolute_value,slowest_without_splitting,fastest_without_splitting,fastest_in_absolute_value_without_splitting)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#ifdef SW
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
#endif
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
#ifdef IMEX
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: HyperbolicityLoss
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)           :: Prim(1:nVar)
REAL,INTENT(OUT),OPTIONAL :: slowest
REAL,INTENT(OUT),OPTIONAL :: fastest
REAL,INTENT(OUT),OPTIONAL :: fastest_in_absolute_value
REAL,INTENT(OUT),OPTIONAL :: slowest_without_splitting
REAL,INTENT(OUT),OPTIONAL :: fastest_without_splitting
REAL,INTENT(OUT),OPTIONAL :: fastest_in_absolute_value_without_splitting
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL                      :: ro, vx, vy, p, ss, ss_without_splitting
#ifdef MULTIFLUID
REAL                      :: roPhi
REAL                      :: pinf
#endif
#ifdef IMEX
REAL             :: problem_if_negative
#endif
!-------------------------------------------------------------------------------!

ro = Prim(1)
#ifdef MOMENTUMINPRIMITIVEVARIABLES
vx = Prim(2)/Prim(1)
vy = Prim(3)/Prim(1)
#else
vx = Prim(2)
vy = Prim(3)
#endif
p  = Prim(4)
#ifdef SW
p  = Kappa*ro**Gmm
#endif
#ifdef MULTIFLUID
roPhi = Prim(5)
Gmm=Get_Gamma(roPhi)
pinf=Get_Pressure_Infinity(roPhi)
#endif


!*Defining the sound speed

#if defined(IMEX)

  IF (( max_ro .LT. ro  ) .OR.  ( min_p .GT. p  )) THEN
    problem_if_negative=0.0
    HyperbolicityLoss=.TRUE.
  ELSE
    problem_if_negative=(max_ro-ro)*(p-min_p)
  END IF
  ss  = SQRT(Gmm*problem_if_negative/(ro*max_ro))

#elif defined(IMEXMOMENTUM)

  ss=0.0

#else

#ifdef MULTIFLUID
  ss  = SQRT(ABS(Gmm*(p+pinf)/ro))
#else
  ss  = SQRT(ABS(Gmm*p/ro))
#endif

#endif

#ifdef RELAXATION
    ss  = ss/EPS_LM
#endif



IF(PRESENT(slowest)) THEN
slowest = vx - ss
END IF

IF(PRESENT(fastest)) THEN
fastest= vx + ss
END IF

IF(PRESENT(fastest_in_absolute_value)) THEN
fastest_in_absolute_value = ABS(vx) + ss
END IF


#ifdef MULTIFLUID
  ss_without_splitting  = SQRT(ABS(Gmm*(p+pinf)/ro))
#else
  ss_without_splitting  = SQRT(ABS(Gmm*p/ro))
#endif

#ifdef RELAXATION
    ss_without_splitting  = ss_without_splitting/EPS_LM
#endif


IF(PRESENT(slowest_without_splitting)) THEN
slowest_without_splitting = vx - ss_without_splitting
END IF

IF(PRESENT(fastest_without_splitting)) THEN
fastest_without_splitting= vx + ss_without_splitting
END IF

IF(PRESENT(fastest_in_absolute_value_without_splitting)) THEN
fastest_in_absolute_value_without_splitting = ABS(vx) + ss_without_splitting
END IF

!-------------------------------------------------------------------------------!
END SUBROUTINE WaveSpeeds1D
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WaveSpeeds2D(Prim,slowestx,fastestx,fastest_in_absolute_valuex,slowesty,fastesty,fastest_in_absolute_valuey &
                          & ,slowestx_without_splitting,fastestx_without_splitting,fastest_in_absolute_valuex_without_splitting,slowesty_without_splitting,fastesty_without_splitting,fastest_in_absolute_valuey_without_splitting)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#ifdef SW
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
#endif
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
#ifdef IMEX
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: HyperbolicityLoss
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT),OPTIONAL :: slowestx
REAL,INTENT(OUT),OPTIONAL :: fastestx
REAL,INTENT(OUT),OPTIONAL :: fastest_in_absolute_valuex
REAL,INTENT(OUT),OPTIONAL :: slowesty
REAL,INTENT(OUT),OPTIONAL :: fastesty
REAL,INTENT(OUT),OPTIONAL :: fastest_in_absolute_valuey
REAL,INTENT(OUT),OPTIONAL :: slowestx_without_splitting
REAL,INTENT(OUT),OPTIONAL :: fastestx_without_splitting
REAL,INTENT(OUT),OPTIONAL :: fastest_in_absolute_valuex_without_splitting
REAL,INTENT(OUT),OPTIONAL :: slowesty_without_splitting
REAL,INTENT(OUT),OPTIONAL :: fastesty_without_splitting
REAL,INTENT(OUT),OPTIONAL :: fastest_in_absolute_valuey_without_splitting
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: ro, vx, vy, p, ss, ss_without_splitting
#ifdef MULTIFLUID
REAL             :: roPhi
REAL             :: pinf
#endif
#ifdef IMEX
REAL             :: problem_if_negative
#endif
!-------------------------------------------------------------------------------!

ro = Prim(1)
#ifdef MOMENTUMINPRIMITIVEVARIABLES
vx = Prim(2)/Prim(1)
vy = Prim(3)/Prim(1)
#else
vx = Prim(2)
vy = Prim(3)
#endif
p  = Prim(4)
#ifdef SW
p  = Kappa*ro**Gmm
#endif
#ifdef MULTIFLUID
roPhi  = Prim(5)
Gmm=Get_Gamma(roPhi)
pinf=Get_Pressure_Infinity(roPhi)
#endif

!*Defining the sound speed

#if defined(IMEX)

  IF (( max_ro .LT. ro  ) .OR.  ( min_p .GT. p  )) THEN
    problem_if_negative=0.0
    HyperbolicityLoss=.TRUE.
  ELSE
    problem_if_negative=(max_ro-ro)*(p-min_p)
  END IF
  ss  = SQRT(Gmm*problem_if_negative/(ro*max_ro))

#elif defined(IMEXMOMENTUM)

  ss=0.0

#else

#ifdef MULTIFLUID
  ss  = SQRT(ABS(Gmm*(p+pinf)/ro))
#else
  ss  = SQRT(ABS(Gmm*p/ro))
#endif

#endif

#ifdef RELAXATION
    ss  = ss/EPS_LM
#endif



IF(PRESENT(slowestx)) THEN
slowestx = vx - ss
END IF

IF(PRESENT(fastestx)) THEN
fastestx= vx + ss
END IF

IF(PRESENT(fastest_in_absolute_valuex)) THEN
fastest_in_absolute_valuex = ABS(vx) + ss
END IF


IF(PRESENT(slowesty)) THEN
slowesty = vy - ss
END IF

IF(PRESENT(fastesty)) THEN
fastesty= vy + ss
END IF

IF(PRESENT(fastest_in_absolute_valuey)) THEN
fastest_in_absolute_valuey = ABS(vy) + ss
END IF


#ifdef MULTIFLUID
  ss_without_splitting  = SQRT(ABS(Gmm*(p+pinf)/ro))
#else
  ss_without_splitting  = SQRT(ABS(Gmm*p/ro))
#endif

#ifdef RELAXATION
    ss_without_splitting  = ss_without_splitting/EPS_LM
#endif


IF(PRESENT(slowestx_without_splitting)) THEN
slowestx_without_splitting = vx - ss_without_splitting
END IF

IF(PRESENT(fastestx_without_splitting)) THEN
fastestx_without_splitting= vx + ss_without_splitting
END IF

IF(PRESENT(fastest_in_absolute_valuex_without_splitting)) THEN
fastest_in_absolute_valuex_without_splitting = ABS(vx) + ss_without_splitting
END IF


IF(PRESENT(slowesty_without_splitting)) THEN
slowesty_without_splitting = vy - ss_without_splitting
END IF

IF(PRESENT(fastesty_without_splitting)) THEN
fastesty_without_splitting= vy + ss_without_splitting
END IF

IF(PRESENT(fastest_in_absolute_valuey_without_splitting)) THEN
fastest_in_absolute_valuey_without_splitting = ABS(vy) + ss_without_splitting
END IF

!-------------------------------------------------------------------------------!
END SUBROUTINE WaveSpeeds2D
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ConsToPrim(Cons, Prim)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DENSITY, MIN_SPEED, MIN_PRESSURE
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Cons(1:nVar)
REAL,INTENT(OUT) :: Prim(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: ro, rovx, rovy, Energy, rot, vx, vy
#ifdef MULTIFLUID
REAL             :: roPhi
REAL             :: pinf
#endif
!-------------------------------------------------------------------------------!

ro     = Cons(1)
rovx   = Cons(2)
rovy   = Cons(3)
Energy = Cons(4)



#ifdef MULTIFLUID
roPhi=Cons(5)
Gmm=Get_Gamma(roPhi)
pinf=Get_Pressure_Infinity(roPhi)
#endif

#ifdef PATANKAR
rot = ro + MIN_DENSITY/ro
#else
#ifdef POSITIVITYCHECKS
IF (ro .LT. MIN_DENSITY) THEN
  ro = MIN_DENSITY
  rovx = 0.
  rovy = 0.
  vx=0.0
  vy=0.0
ELSE
  vx=rovx/ro
  vy=rovy/ro
END IF
#else

vx=rovx/ro
vy=rovy/ro

#endif
rot = ro
#endif


Prim(1) = ro
#ifdef MOMENTUMINPRIMITIVEVARIABLES
Prim(2) = rovx
Prim(3) = rovy
#else
Prim(2) = rovx/rot
Prim(3) = rovy/rot
#endif

#ifdef RELAXATION
Prim(4) = (Gmm-1.0)*( Energy-EPS_LM**2*0.5*ro*(vx**2+vy**2) )
#else
Prim(4) = (Gmm-1.0)*( Energy-0.5*ro*(vx**2+vy**2) )
#endif


#ifdef MULTIFLUID
Prim(4) = Prim(4)-Gmm*pinf
#endif


#ifdef POSITIVITYCHECKS
IF (Prim(4) .LT. MIN_PRESSURE) THEN
  Prim(4) = MIN_PRESSURE
  Prim(2) = 0.
  Prim(3) = 0.
END IF
#endif

#ifdef MULTIFLUID

#ifdef PRIMITIVEFORMULATIONLEVELSET
      Prim(5)=roPhi/ro !*Primitive is only phi in this case
#else
      Prim(5)=roPhi
#endif

#endif




!-------------------------------------------------------------------------------!
END SUBROUTINE ConsToPrim
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE PrimToCons(Prim, Cons)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DENSITY, MIN_SPEED, MIN_PRESSURE
#ifdef SW
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
#endif
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: Cons(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: ro, vx, vy, p
#ifdef MULTIFLUID
REAL             :: roPhi
REAL             :: pinf
#endif
!-------------------------------------------------------------------------------!

ro  = Prim(1)
#ifdef MOMENTUMINPRIMITIVEVARIABLES
vx  = Prim(2)/Prim(1)
vy  = Prim(3)/Prim(1)
#else
vx  = Prim(2)
vy  = Prim(3)
#endif
p   = Prim(4)
#ifdef SW
p   = Kappa*ro**Gmm
#endif
#ifdef MULTIFLUID
roPhi=Prim(5)
Gmm =Get_Gamma(roPhi)
pinf=Get_Pressure_Infinity(roPhi)
#endif


#ifdef PATANKAR

#else
#ifdef POSITIVITYCHECKS
IF (ro .LT. MIN_DENSITY) THEN
 ro = MIN_DENSITY
 vx=0.
 vy=0.
END IF
#endif
#endif

#ifdef POSITIVITYCHECKS
IF (p .LT. MIN_PRESSURE) THEN
  p = MIN_PRESSURE
  vx = 0.
  vy = 0.
END IF
#endif

Cons(1) = ro
Cons(2) = ro*vx
Cons(3) = ro*vy
#ifdef RELAXATION
Cons(4) = p/(Gmm-1.0)+EPS_LM**2*0.5*ro*(vx**2+vy**2)
#else
Cons(4) = p/(Gmm-1.0)+0.5*ro*(vx**2+vy**2)
#endif

#ifdef MULTIFLUID
Cons(4) = Cons(4) + Gmm/(Gmm-1.0)*pinf
#endif


#ifdef MULTIFLUID


#ifdef PRIMITIVEFORMULATIONLEVELSET
      Cons(5)=roPhi*ro !*Primitive is phi in this case
#else
      Cons(5)=roPhi
#endif

#endif


!-------------------------------------------------------------------------------!
END SUBROUTINE PrimToCons
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE EvaluateFlux1D(Prim,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DENSITY, MIN_SPEED, MIN_PRESSURE
#ifdef SW
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
#endif
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: Flux(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: ro, vx, vy, p, Energy
#ifdef MULTIFLUID
REAL             :: roPhi
REAL             :: pinf
#endif
!-------------------------------------------------------------------------------!

ro = Prim(1)
#ifdef MOMENTUMINPRIMITIVEVARIABLES
vx = Prim(2)/Prim(1)
vy = Prim(3)/Prim(1)
#else
vx = Prim(2)
vy = Prim(3)
#endif
p  = Prim(4)
#ifdef SW
p  = Kappa*ro**Gmm
#endif
#ifdef MULTIFLUID
roPhi=Prim(5)
Gmm=Get_Gamma(roPhi)
pinf=Get_Pressure_Infinity(roPhi)
#endif

#ifdef PATANKAR

#else
#ifdef POSITIVITYCHECKS
IF (ro .LT. MIN_DENSITY) THEN
 ro = MIN_DENSITY
 vx=0.
 vy=0.
END IF
#endif
#endif

#ifdef POSITIVITYCHECKS
IF (p .LT. MIN_PRESSURE) THEN
  p = MIN_PRESSURE
  vx = 0.
  vy = 0.
END IF
#endif

#ifdef RELAXATION
Energy = p/(Gmm-1.0)+EPS_LM**2*0.5*ro*(vx**2+vy**2)
#else
Energy = p/(Gmm-1.0)+0.5*ro*(vx**2+vy**2)
#endif

#ifdef MULTIFLUID
Energy = Energy+Gmm/(Gmm-1.0)*pinf
#endif

Flux(1) = ro*vx
#ifdef RELAXATION
Flux(2) = ro*vx**2 + p/EPS_LM**2
#else
Flux(2) = ro*vx**2 + p
#endif
Flux(3) = ro*vx*vy
Flux(4) = vx*(Energy+p)
!*NB: Reference:
!*Eleuterio F. Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics - A Practical Introduction
!*3.2.4 The Split Three–Dimensional Riemann Problem

#ifdef MULTIFLUID
#ifdef LEVELSETINU

#ifdef PRIMITIVEFORMULATIONLEVELSET
Flux(5)=vx*roPhi*ro !*Primitive is phi in this case
#else
Flux(5)=vx*roPhi
#endif

#else
Flux(5)=0.0
#endif
#endif
!-------------------------------------------------------------------------------!
END SUBROUTINE EvaluateFlux1D
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE RiemannSolver(PrimL,PrimR,NormVect,TangVect,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: WhichRiemannSolver
USE exact_riemann_mod,      ONLY: exact_riemann
#ifdef SW
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
#endif
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: PrimL(1:nVar,1:nGPs)
REAL,INTENT(IN)  :: PrimR(1:nVar,1:nGPs)
REAL,INTENT(IN)  :: NormVect(1:nDims)
REAL,INTENT(IN)  :: TangVect(1:nDims)
REAL,INTENT(OUT) :: Flux(1:nVar,1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: PrimLL(1:nVar,1:nGPs), PrimRR(1:nVar,1:nGPs)
REAL             :: ConsLL(1:nVar,1:nGPs), ConsRR(1:nVar,1:nGPs)
INTEGER          :: iGP
!-------------------------------------------------------------------------------!
! >> LOCAL FOR EXACT RIEMANN SOLVER                                             !
!-------------------------------------------------------------------------------!
REAL                        :: PrimStar(1:nVar)
REAL, PARAMETER             :: s=0.0
REAL                        :: al, ar
REAL                        :: pl, pr
REAL                        :: rho_star_l,rho_star_r
REAL                        :: speedl, speedr
REAL                        :: pm, um
REAL                        :: u_norm_l, u_norm_r 
REAL                        :: u_tan_l,  u_tan_r
REAL, DIMENSION(3+nDims)    :: vstar
REAL, DIMENSION(2+nDims)    :: w
!-------------------------------------------------------------------------------!



DO iGP=1,nGPs
  ! Rotating the vector quantities       !
  PrimLL(1,iGP) = PrimL(1,iGP)
  PrimLL(2,iGP) = NormVect(1)*PrimL(2,iGP) + NormVect(2)*PrimL(3,iGP)
  PrimLL(3,iGP) = TangVect(1)*PrimL(2,iGP) + TangVect(2)*PrimL(3,iGP)
  PrimLL(4,iGP) = PrimL(4,iGP)
#ifdef MULTIFLUID
  PrimLL(5,iGP) = PrimL(5,iGP)
#endif

  PrimRR(1,iGP) = PrimR(1,iGP)
  PrimRR(2,iGP) = NormVect(1)*PrimR(2,iGP) + NormVect(2)*PrimR(3,iGP)
  PrimRR(3,iGP) = TangVect(1)*PrimR(2,iGP) + TangVect(2)*PrimR(3,iGP)
  PrimRR(4,iGP) = PrimR(4,iGP)
#ifdef MULTIFLUID
  PrimRR(5,iGP) = PrimR(5,iGP)
#endif

  CALL PrimToCons(PrimLL(1:nVar,iGP),ConsLL(1:nVar,iGP))
  CALL PrimToCons(PrimRR(1:nVar,iGP),ConsRR(1:nVar,iGP))

  SELECT CASE(WhichRiemannSolver)
    CASE(0) !*Rusanov
      CALL RiemannSolverByRusanov(&
        ConsLL(1:nVar,iGP),ConsRR(1:nVar,iGP),&
        PrimLL(1:nVar,iGP),PrimRR(1:nVar,iGP),Flux(1:nVar,iGP))
    CASE(-1)
      CALL CentralUpwind(&
        ConsLL(1:nVar,iGP),ConsRR(1:nVar,iGP),&
        PrimLL(1:nVar,iGP),PrimRR(1:nVar,iGP),Flux(1:nVar,iGP))
    CASE(1)
      CALL LowDissipationCentralUpwind(&
        ConsLL(1:nVar,iGP),ConsRR(1:nVar,iGP),&
        PrimLL(1:nVar,iGP),PrimRR(1:nVar,iGP),Flux(1:nVar,iGP))
    CASE(2)
      CALL LowDissipationCentralUpwindNewAntidiffusionTerm(&
        ConsLL(1:nVar,iGP),ConsRR(1:nVar,iGP),&
        PrimLL(1:nVar,iGP),PrimRR(1:nVar,iGP),Flux(1:nVar,iGP))
    CASE(3) !*Exact

#ifdef RELAXATION
      IF (EPS_LM .NE. 1.0) THEN
        PRINT*, "Exact Riemann solver needs to be adjusted for relaxation"
        STOP
      END IF
#endif

#ifdef MULTIFLUID
      PRINT*, "You shouldn't be here in RiemannSolver"
      STOP
#endif

      CALL exact_riemann(Gmm,PrimLL(1,iGP),PrimLL(2,iGP),PrimLL(4,iGP),PrimRR(1,iGP),PrimRR(2,iGP),PrimRR(4,iGP),PrimStar(1),PrimStar(2),PrimStar(4))

      IF (PrimStar(2) .GT. 0.0) THEN
          PrimStar(3)=PrimLL(3,iGP)
      ELSE
          PrimStar(3)=PrimRR(3,iGP)
      ENDIF


      CALL EvaluateFlux1D(PrimStar(1:nVar),Flux(1:nVar,iGP))

    CASE DEFAULT
      PRINT*, "Riemann Solver not defined"
      PRINT*, "Riemann Solver was", WhichRiemannSolver
      STOP
  END SELECT

  ! Rotating back the momentum components
  Flux(2:3,iGP) = NormVect(1:nDims)*Flux(2,iGP) &
                + TangVect(1:nDims)*Flux(3,iGP)
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE RiemannSolver
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE RiemannSolverByRusanov(ConsL,ConsR,PrimL,PrimR,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: ConsL(1:nVar), ConsR(1:nVar)
REAL,INTENT(IN)  :: PrimL(1:nVar), PrimR(1:nVar)
REAL,INTENT(OUT) :: Flux(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: FluxL(1:nVar), FluxR(1:nVar)
REAL             :: LambdaMax, fastest_in_absolute_valueL, fastest_in_absolute_valueR
!-------------------------------------------------------------------------------!

CALL EvaluateFlux1D(PrimL,FluxL)
CALL EvaluateFlux1D(PrimR,FluxR)
CALL WaveSpeeds1D(PrimL,fastest_in_absolute_value_without_splitting=fastest_in_absolute_valueL)
CALL WaveSpeeds1D(PrimR,fastest_in_absolute_value_without_splitting=fastest_in_absolute_valueR)

LambdaMax = MAX(fastest_in_absolute_valueL,fastest_in_absolute_valueR)

Flux = 0.5*((FluxL + FluxR) - LambdaMax*(ConsR - ConsL))

!-------------------------------------------------------------------------------!
END SUBROUTINE RiemannSolverByRusanov
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE CentralUpwind(ConsL,ConsR,PrimL,PrimR,Flux)
!-------------------------------------------------------------------------------!
!*New Low-Dissipation Central-Upwind Schemes, Alexander Kurganov, Ruixiao Xin
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: EPS_AVOIDING_DIVISION_BY_ZERO_IN_CU
USE MOD_Reconstruction,     ONLY: MINMOD
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: ConsL(1:nVar), ConsR(1:nVar)
REAL,INTENT(IN)  :: PrimL(1:nVar), PrimR(1:nVar)
REAL,INTENT(OUT) :: Flux(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: FluxL(1:nVar), FluxR(1:nVar)
REAL             :: slowestL, fastestL, slowestR, fastestR
REAL             :: SL, SM, SR
REAL             :: dC(1:nVar)                         !*Built-in antidiffusion term
REAL             :: Cstar(1:nVar)                      !*Intermediate value of C
INTEGER          :: iVar
!-------------------------------------------------------------------------------!

CALL EvaluateFlux1D(PrimL,FluxL)
CALL EvaluateFlux1D(PrimR,FluxR)

CALL WaveSpeeds1D(PrimL,slowest_without_splitting=slowestL,fastest_without_splitting=fastestL)
CALL WaveSpeeds1D(PrimR,slowest_without_splitting=slowestR,fastest_without_splitting=fastestR)
SL=MIN(slowestL,slowestR,-EPS_AVOIDING_DIVISION_BY_ZERO_IN_CU)
SR=MAX(fastestL,fastestR,+EPS_AVOIDING_DIVISION_BY_ZERO_IN_CU)

Flux = (SR * FluxL - SL * FluxR + SL * SR * (ConsR - ConsL)) / (SR - SL)

!-------------------------------------------------------------------------------!
END SUBROUTINE CentralUpwind
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE LowDissipationCentralUpwind(ConsL,ConsR,PrimL,PrimR,Flux)
!-------------------------------------------------------------------------------!
!*New Low-Dissipation Central-Upwind Schemes, Alexander Kurganov, Ruixiao Xin
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: EPS_AVOIDING_DIVISION_BY_ZERO_IN_CU
USE MOD_Reconstruction,     ONLY: MINMOD
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: ConsL(1:nVar), ConsR(1:nVar)
REAL,INTENT(IN)  :: PrimL(1:nVar), PrimR(1:nVar)
REAL,INTENT(OUT) :: Flux(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: FluxL(1:nVar), FluxR(1:nVar)
REAL             :: slowestL, fastestL, slowestR, fastestR
REAL             :: SL, SM, SR
REAL             :: dC(1:nVar)                         !*Built-in antidiffusion term
REAL             :: Cstar(1:nVar)                      !*Intermediate value of C
INTEGER          :: iVar
!-------------------------------------------------------------------------------!

CALL EvaluateFlux1D(PrimL,FluxL)
CALL EvaluateFlux1D(PrimR,FluxR)

CALL WaveSpeeds1D(PrimL,slowest_without_splitting=slowestL,fastest_without_splitting=fastestL)
CALL WaveSpeeds1D(PrimR,slowest_without_splitting=slowestR,fastest_without_splitting=fastestR)
SL=MIN(slowestL,slowestR,-EPS_AVOIDING_DIVISION_BY_ZERO_IN_CU)
SR=MAX(fastestL,fastestR,+EPS_AVOIDING_DIVISION_BY_ZERO_IN_CU)

!*Built-in antidiffusion term
Cstar = ( SR*ConsR - SL*ConsL - ( FluxR - FluxL ) ) / ( SR - SL )
dC = 0.0
Do iVar=1,nVar
  dC(iVar)=MINMOD( Cstar(iVar)-ConsL(iVar), ConsR(iVar)-Cstar(iVar) )
END DO

Flux = (SR * FluxL - SL * FluxR + SL * SR * (ConsR - ConsL - dC)) / (SR - SL)

!-------------------------------------------------------------------------------!
END SUBROUTINE LowDissipationCentralUpwind
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE LowDissipationCentralUpwindNewAntidiffusionTerm(ConsL,ConsR,PrimL,PrimR,Flux)
!-------------------------------------------------------------------------------!
!*New Low-Dissipation Central-Upwind Schemes, Alexander Kurganov, Ruixiao Xin
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: EPS_AVOIDING_DIVISION_BY_ZERO_IN_CU
USE MOD_Reconstruction,     ONLY: MINMOD
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: ConsL(1:nVar), ConsR(1:nVar)
REAL,INTENT(IN)  :: PrimL(1:nVar), PrimR(1:nVar)
REAL,INTENT(OUT) :: Flux(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: FluxL(1:nVar), FluxR(1:nVar)
REAL             :: slowestL, fastestL, slowestR, fastestR
REAL             :: SL, SM, SR
REAL             :: dC(1:nVar)                         !*Built-in antidiffusion term
REAL             :: Cstar(1:nVar)                      !*Intermediate value of C
INTEGER          :: iVar
REAL             :: uxstar, uystar, alphastar, w1, w2, Q(1:nVar), qro, ap, am, www, dE, wl, wr, drho, ustar, dmystar, dmy, dmx, dmxstar, rhostar
!-------------------------------------------------------------------------------!

CALL EvaluateFlux1D(PrimL,FluxL)
CALL EvaluateFlux1D(PrimR,FluxR)

CALL WaveSpeeds1D(PrimL,slowest_without_splitting=slowestL,fastest_without_splitting=fastestL)
CALL WaveSpeeds1D(PrimR,slowest_without_splitting=slowestR,fastest_without_splitting=fastestR)
SL=MIN(slowestL,slowestR,-EPS_AVOIDING_DIVISION_BY_ZERO_IN_CU)
SR=MAX(fastestL,fastestR,+EPS_AVOIDING_DIVISION_BY_ZERO_IN_CU)
am=SL
ap=SR
www=SR-SL

          rhostar=((SR*ConsR(1)-SL*ConsL(1))-(FluxR(1)-FluxL(1)))/www
          dmxstar =((SR*ConsR(2)-SL*ConsL(2))-(FluxR(2)-FluxL(2)))/www
          dmystar =((SR*ConsR(3)-SL*ConsL(3))-(FluxR(3)-FluxL(3)))/www
          ustar=dmxstar/rhostar

          wl=-(SL-ustar)*(ConsL(1)-rhostar)
          wr= (SR-ustar)*(rhostar-ConsR(1))
          drho=MINMOD(wl,wr)
      
          wl=-(SL-ustar)*(ConsL(3)-dmystar)
          wr= (SR-ustar)*(dmystar-ConsR(3))
          dmy=MINMOD(wl,wr)

          dE=0.50*((dmystar-dmy/(SR-ustar))**2/(rhostar-drho/(SR-ustar))-(dmystar-dmy/(SL-ustar))**2/(rhostar-drho/(SL-ustar)))/(SR-SL)

          if(ustar>0.0)then
           drho=drho*SL/(SL-ustar)
           dmy =dmy*SL/(SL-ustar)
           dE  =dE*SL*(SR-ustar)
          else
           drho=drho*SR/(SR-ustar)
           dmy =dmy*SR/(SR-ustar)
           dE  =dE*SR*(SL-ustar)    
          endif

          dmx=ustar*drho
          dE=dE+0.50*ustar**2*drho

          Q(1) = drho
          Q(2) = dmx
          Q(3) = dmy
          Q(4) = dE

          Flux = ( ap*FluxL - am*FluxR +ap*am*(ConsR-ConsL) ) / (ap-am) - Q





#ifdef OLDEVOLVERHOWITHSAMESCHEME
!*Built-in antidiffusion term
Cstar = ( SR*ConsR - SL*ConsL - ( FluxR - FluxL ) ) / ( SR - SL )

dC = 0.0
iVar = 1
dC(iVar)=MINMOD( Cstar(iVar)-ConsL(iVar), ConsR(iVar)-Cstar(iVar) )
Flux(iVar) = (SR * FluxL(iVar) - SL * FluxR(iVar) + SL * SR * (ConsR(iVar) - ConsL(iVar) - dC(iVar))) / (SR - SL)
#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE LowDissipationCentralUpwindNewAntidiffusionTerm
!===============================================================================!
!
!
!
!===============================================================================!
#ifndef MULTIFLUID
SUBROUTINE LeftEigenvectors(Cons,NormVect,LMat)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: nDims
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: Cons(1:nVar)
REAL,INTENT(IN)    :: NormVect(1:nDims)
REAL,INTENT(INOUT) :: LMat(1:nVar,1:nVar)
!-------------------------------------------------------------------------------!
REAL               :: v_nn(nDims)

REAL               :: velocity(nDims)
REAL               :: TangVect(nDims) 

REAL               :: ek, v_n, v_t, ht, speed_of_sound, pi_rho, pi_E, xi
REAL               :: ro, E, p, eps, vx, vy

REAL               :: RMat(1:nVar,1:nVar)

#if defined(RELAXATION)
! PRINT*, "Stopping in LeftEigenvectors. I do not know what to do in IMEX case"
#else

velocity  = Cons(2:3)/Cons(1)                                       ! > u vector

v_n = SUM(velocity*NormVect)                                                ! > u.n

ek = 0.50 * SUM(velocity**2)                                         ! > Kinetic Energy

!EOS
ro             = Cons(1)
vx             = Cons(2)/Cons(1)
vy             = Cons(3)/Cons(1)
eps            = Cons(4)/Cons(1) - 0.50*( SUM(Cons(2:3)**2) )/Cons(1)**2
p              = (Gmm-1.0)*( Cons(4)-0.5*ro*(vx**2+vy**2) )
speed_of_sound = SQRT(ABS(Gmm*p/ro)) !*ALERT: Absolute value
E              = Cons(4)
ht             = (E+p)/ro                                                ! > Total Enthalpy

TangVect       = 0.
TangVect(1)    = NormVect(2); TangVect(2) = -NormVect(1)                               ! > nn_p.nn = 0
v_t            = SUM(velocity*TangVect)                                              ! > u_p = u.nn_p

pi_rho         = (Gmm-1.0)*0.5*(SUM(Cons(2:3)**2))/Cons(1)**2                     ! > Pressure Derivative // rho
pi_E           = Gmm-1.                                                           ! > Pressure Derivative // rho*E
xi             = 2.0*ek - pi_rho/pi_E

LMat = 0.0

LMat(1, 1)         =  speed_of_sound*v_n + pi_rho
LMat(1, 2:nDims+1) = -speed_of_sound*NormVect  - pi_E*velocity
LMat(1, 4)         =  pi_E

LMat(2, 1)         = -2.0 * (pi_rho - speed_of_sound**2)
LMat(2, 2:nDims+1) =  2.0 *  pi_E*velocity
LMat(2, 4)         = -2.0 *  pi_E

LMat(3, 1)         = -speed_of_sound*v_n + pi_rho
LMat(3, 2:nDims+1) =  speed_of_sound*NormVect  - pi_E*velocity
LMat(3, 4)         =  pi_E

LMat(4, 1)         = -2.0*speed_of_sound**2 * v_t
LMat(4, 2:nDims+1) =  2.0*speed_of_sound**2 * TangVect
LMat(4, 4)         =  0.0

LMat = LMat / (2.0 * (speed_of_sound**2) )


! CALL RightEigenvectors(Cons,NormVect,RMat)
! CALL M44INV(RMat,LMat)


#ifdef MULTIFLUID
! PRINT*, "I do not know the LeftEigenvectors for multifluid"
! STOP
#endif

#endif




CONTAINS

  SUBROUTINE M44INV (A, AINV)
    IMPLICIT NONE
    REAL, DIMENSION(4,4), INTENT(IN)  :: A
    REAL, DIMENSION(4,4), INTENT(OUT) :: AINV

    REAL :: DET
    REAL, DIMENSION(4,4) :: COFACTOR

    ! Compute the determinant of the 4x4 matrix
    DET =   A(1,1)*(A(2,2)*A(3,3)*A(4,4) + A(2,3)*A(3,4)*A(4,2) + A(2,4)*A(3,2)*A(4,3) &
              - A(2,4)*A(3,3)*A(4,2) - A(2,2)*A(3,4)*A(4,3) - A(2,3)*A(3,2)*A(4,4)) &
          - A(1,2)*(A(2,1)*A(3,3)*A(4,4) + A(2,3)*A(3,4)*A(4,1) + A(2,4)*A(3,1)*A(4,3) &
              - A(2,4)*A(3,3)*A(4,1) - A(2,1)*A(3,4)*A(4,3) - A(2,3)*A(3,1)*A(4,4)) &
          + A(1,3)*(A(2,1)*A(3,2)*A(4,4) + A(2,2)*A(3,4)*A(4,1) + A(2,4)*A(3,1)*A(4,2) &
              - A(2,4)*A(3,2)*A(4,1) - A(2,1)*A(3,4)*A(4,2) - A(2,2)*A(3,1)*A(4,4)) &
          - A(1,4)*(A(2,1)*A(3,2)*A(4,3) + A(2,2)*A(3,3)*A(4,1) + A(2,3)*A(3,1)*A(4,2) &
              - A(2,3)*A(3,2)*A(4,1) - A(2,1)*A(3,3)*A(4,2) - A(2,2)*A(3,1)*A(4,3))

    IF (DET == 0.0) THEN
        PRINT *, "The matrix is singular and cannot be inverted."
        RETURN
    END IF

    ! Compute the cofactors of the 4x4 matrix
    COFACTOR(1,1) = +(A(2,2)*A(3,3)*A(4,4) + A(2,3)*A(3,4)*A(4,2) + A(2,4)*A(3,2)*A(4,3) &
                   - A(2,4)*A(3,3)*A(4,2) - A(2,2)*A(3,4)*A(4,3) - A(2,3)*A(3,2)*A(4,4))
    COFACTOR(1,2) = -(A(2,1)*A(3,3)*A(4,4) + A(2,3)*A(3,4)*A(4,1) + A(2,4)*A(3,1)*A(4,3) &
                   - A(2,4)*A(3,3)*A(4,1) - A(2,1)*A(3,4)*A(4,3) - A(2,3)*A(3,1)*A(4,4))
    COFACTOR(1,3) = +(A(2,1)*A(3,2)*A(4,4) + A(2,2)*A(3,4)*A(4,1) + A(2,4)*A(3,1)*A(4,2) &
                   - A(2,4)*A(3,2)*A(4,1) - A(2,1)*A(3,4)*A(4,2) - A(2,2)*A(3,1)*A(4,4))
    COFACTOR(1,4) = -(A(2,1)*A(3,2)*A(4,3) + A(2,2)*A(3,3)*A(4,1) + A(2,3)*A(3,1)*A(4,2) &
                   - A(2,3)*A(3,2)*A(4,1) - A(2,1)*A(3,3)*A(4,2) - A(2,2)*A(3,1)*A(4,3))

    COFACTOR(2,1) = -(A(1,2)*A(3,3)*A(4,4) + A(1,3)*A(3,4)*A(4,2) + A(1,4)*A(3,2)*A(4,3) &
                   - A(1,4)*A(3,3)*A(4,2) - A(1,2)*A(3,4)*A(4,3) - A(1,3)*A(3,2)*A(4,4))
    COFACTOR(2,2) = +(A(1,1)*A(3,3)*A(4,4) + A(1,3)*A(3,4)*A(4,1) + A(1,4)*A(3,1)*A(4,3) &
                   - A(1,4)*A(3,3)*A(4,1) - A(1,1)*A(3,4)*A(4,3) - A(1,3)*A(3,1)*A(4,4))
    COFACTOR(2,3) = -(A(1,1)*A(3,2)*A(4,4) + A(1,2)*A(3,4)*A(4,1) + A(1,4)*A(3,1)*A(4,2) &
                   - A(1,4)*A(3,2)*A(4,1) - A(1,1)*A(3,4)*A(4,2) - A(1,2)*A(3,1)*A(4,4))
    COFACTOR(2,4) = +(A(1,1)*A(3,2)*A(4,3) + A(1,2)*A(3,3)*A(4,1) + A(1,3)*A(3,1)*A(4,2) &
                   - A(1,3)*A(3,2)*A(4,1) - A(1,1)*A(3,3)*A(4,2) - A(1,2)*A(3,1)*A(4,3))

    COFACTOR(3,1) = +(A(1,2)*A(2,3)*A(4,4) + A(1,3)*A(2,4)*A(4,2) + A(1,4)*A(2,2)*A(4,3) &
                   - A(1,4)*A(2,3)*A(4,2) - A(1,2)*A(2,4)*A(4,3) - A(1,3)*A(2,2)*A(4,4))
    COFACTOR(3,2) = -(A(1,1)*A(2,3)*A(4,4) + A(1,3)*A(2,4)*A(4,1) + A(1,4)*A(2,1)*A(4,3) &
                   - A(1,4)*A(2,3)*A(4,1) - A(1,1)*A(2,4)*A(4,3) - A(1,3)*A(2,1)*A(4,4))
    COFACTOR(3,3) = +(A(1,1)*A(2,2)*A(4,4) + A(1,2)*A(2,4)*A(4,1) + A(1,4)*A(2,1)*A(4,2) &
                   - A(1,4)*A(2,2)*A(4,1) - A(1,1)*A(2,4)*A(4,2) - A(1,2)*A(2,1)*A(4,4))
    COFACTOR(3,4) = -(A(1,1)*A(2,2)*A(4,3) + A(1,2)*A(2,3)*A(4,1) + A(1,3)*A(2,1)*A(4,2) &
                   - A(1,3)*A(2,2)*A(4,1) - A(1,1)*A(2,3)*A(4,2) - A(1,2)*A(2,1)*A(4,3))

    COFACTOR(4,1) = -(A(1,2)*A(2,3)*A(3,4) + A(1,3)*A(2,4)*A(3,2) + A(1,4)*A(2,2)*A(3,3) &
                   - A(1,4)*A(2,3)*A(3,2) - A(1,2)*A(2,4)*A(3,3) - A(1,3)*A(2,2)*A(3,4))
    COFACTOR(4,2) = +(A(1,1)*A(2,3)*A(3,4) + A(1,3)*A(2,4)*A(3,1) + A(1,4)*A(2,1)*A(3,3) &
                   - A(1,4)*A(2,3)*A(3,1) - A(1,1)*A(2,4)*A(3,3) - A(1,3)*A(2,1)*A(3,4))
    COFACTOR(4,3) = -(A(1,1)*A(2,2)*A(3,4) + A(1,2)*A(2,4)*A(3,1) + A(1,4)*A(2,1)*A(3,2) &
                   - A(1,4)*A(2,2)*A(3,1) - A(1,1)*A(2,4)*A(3,2) - A(1,2)*A(2,1)*A(3,4))
    COFACTOR(4,4) = +(A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) &
                   - A(1,3)*A(2,2)*A(3,1) - A(1,1)*A(2,3)*A(3,2) - A(1,2)*A(2,1)*A(3,3))

    ! Transpose the cofactor matrix and divide by the determinant to get the inverse
    AINV = TRANSPOSE(COFACTOR) / DET

  END SUBROUTINE M44INV


!-------------------------------------------------------------------------------!
END SUBROUTINE LeftEigenvectors
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifndef MULTIFLUID
SUBROUTINE RightEigenvectors(Cons,NormVect,RMat)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: nDims
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: Cons(1:nVar)
REAL,INTENT(IN)    :: NormVect(nDims)
REAL,INTENT(INOUT) :: RMat(1:nVar,1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: velocity(nDims)
REAL               :: TangVect(nDims)

REAL :: ek, v_n, v_t, ht, speed_of_sound, pi_rho, pi_E, xi
REAL :: ro, E, p, eps, vx, vy

#if defined(RELAXATION)
! PRINT*, "Stopping in RightEigenvectors. I do not know what to do in IMEX case"
#else

velocity  = Cons(2:3)/Cons(1)                                       ! > velocity vector

v_n = SUM(velocity*NormVect)                                        ! > normal velocity

ek = 0.5 * SUM(velocity**2)                                         ! > Kinetic Energy (up to density)

!EOS
ro            = Cons(1)
vx             = Cons(2)/Cons(1)
vy             = Cons(3)/Cons(1)
eps            = Cons(4)/Cons(1) - 0.5*( SUM(Cons(2:3)**2) )/Cons(1)**2
p              = (Gmm-1.0)*( Cons(4)-0.5*ro*(vx**2+vy**2) )
speed_of_sound = SQRT(ABS(Gmm*p/ro)) !*ALERT: Absolute value
E              = Cons(4)
ht             = (E+p)/ro                                                    ! > Total Enthalpy

TangVect       = 0.
TangVect(1)    = NormVect(2); TangVect(2) = -NormVect(1)                               ! > nn_p.nn = 0
v_t            = SUM(velocity*TangVect)                                              ! > u_p = u.nn_p

pi_rho         = (Gmm-1.0)*0.5*(SUM(Cons(2:3)**2))/Cons(1)**2                     ! > Pressure Derivative // rho
pi_E           = Gmm-1.                                                           ! > Pressure Derivative // rho*E
xi             = 2.*ek - pi_rho/pi_E

RMat = 0.0

RMat(1,         1) = 1.
RMat(2:nDims+1, 1) = velocity - speed_of_sound*NormVect
RMat(4,         1) = ht - speed_of_sound*v_n

RMat(1,         2) = 1.
RMat(2:nDims+1, 2) = velocity
RMat(4,         2) = xi

RMat(1,         3) = 1.
RMat(2:nDims+1, 3) = velocity + speed_of_sound*NormVect
RMat(4,         3) = ht + speed_of_sound*v_n

RMat(1,         4) = 0.
RMat(2:nDims+1, 4) = TangVect
RMat(4,         4) = v_t

#if(1==0)
RMat=0.0
RMat(1,1) = 1.0
RMat(2,2) = 1.0
RMat(3,3) = 1.0 
RMat(4,4) = 1.0 
#endif

#ifdef MULTIFLUID
! PRINT*, "I do not know the RightEigenvectors for multifluid"
! STOP
#endif

#endif


!-------------------------------------------------------------------------------!
END SUBROUTINE RightEigenvectors
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifndef MULTIFLUID
SUBROUTINE RightEigenvectorsObtainedFromRotationalInvariance(Cons,NormVect,RMat)
!*------------------------------------------------------------------------------!
!*Toro's book
!*Proposition 3.15 (Rotational Invariance)
!*K(U,theta)=T^{-1}(theta)K(T(theta)U)
!*where K is the matrix of the eigenvectors of the F flux (first component) of the Euler equations and T is the rotation matrix
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: nDims
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: Cons(1:nVar)
REAL,INTENT(IN)    :: NormVect(nDims)
REAL,INTENT(INOUT) :: RMat(1:nVar,1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: TMat(1:nVar,1:nVar)
REAL               :: invTMat(1:nVar,1:nVar)
REAL               :: KMat(1:nVar,1:nVar)
REAL               :: TU(1:nVar)
REAL               :: theta

#if defined(RELAXATION)
! PRINT*, "Stopping in RightEigenvectorsObtainedFromRotationalInvariance. I do not know what to do in IMEX case"
#else

theta = ATAN2(NormVect(2), NormVect(1))

TMat(1:4,1:4)=RESHAPE([1.0,          0.0,          0.0,               0.0, &
          &   0.0,   COS(theta),  -SIN(theta),               0.0, &
          &   0.0,   SIN(theta),   COS(theta),               0.0, &
          &   0.0,          0.0,          0.0,           1.0],[4,4])

invTMat(1:4,1:4)=RESHAPE([1.0,          0.0,          0.0,               0.0, &
          &   0.0,   COS(theta),   SIN(theta),               0.0, &
          &   0.0,  -SIN(theta),   COS(theta),               0.0, &
          &   0.0,          0.0,          0.0,           1.0],[4,4])

TU(1:nVar)=MATMUL(TMat(1:nVar,1:nVar),Cons(1:nVar))

CALL RightEigenvectorsXDirection(TU(1:nVar),KMat(1:nVar,1:nVar))

RMat(1:nVar,1:nVar)=MATMUL(invTMat(1:nVar,1:nVar),KMat(1:nVar,1:nVar))

#endif


!-------------------------------------------------------------------------------!
END SUBROUTINE RightEigenvectorsObtainedFromRotationalInvariance
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifndef MULTIFLUID
SUBROUTINE LeftEigenvectorsObtainedFromRotationalInvariance(Cons,NormVect,LMat)
!*------------------------------------------------------------------------------!
!*Toro's book
!*Proposition 3.15 (Rotational Invariance)
!*K(U,theta)=T^{-1}(theta)K(T(theta)U)
!*where K is the matrix of the eigenvectors of the F flux (first component) of the Euler equations and T is the rotation matrix
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: nDims
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: Cons(1:nVar)
REAL,INTENT(IN)    :: NormVect(1:nDims)
REAL,INTENT(INOUT) :: LMat(1:nVar,1:nVar)
!-------------------------------------------------------------------------------!
REAL               :: RMat(1:nVar,1:nVar)

#if defined(RELAXATION)
! PRINT*, "Stopping in LeftEigenvectorsObtainedFromRotationalInvariance. I do not know what to do in IMEX case"
#else

CALL RightEigenvectorsObtainedFromRotationalInvariance(Cons,NormVect,RMat)
CALL M44INV(RMat,LMat)

#endif


CONTAINS

  SUBROUTINE M44INV (A, AINV)
    IMPLICIT NONE
    REAL, DIMENSION(4,4), INTENT(IN)  :: A
    REAL, DIMENSION(4,4), INTENT(OUT) :: AINV

    REAL :: DET
    REAL, DIMENSION(4,4) :: COFACTOR

    ! Compute the determinant of the 4x4 matrix
    DET =   A(1,1)*(A(2,2)*A(3,3)*A(4,4) + A(2,3)*A(3,4)*A(4,2) + A(2,4)*A(3,2)*A(4,3) &
              - A(2,4)*A(3,3)*A(4,2) - A(2,2)*A(3,4)*A(4,3) - A(2,3)*A(3,2)*A(4,4)) &
          - A(1,2)*(A(2,1)*A(3,3)*A(4,4) + A(2,3)*A(3,4)*A(4,1) + A(2,4)*A(3,1)*A(4,3) &
              - A(2,4)*A(3,3)*A(4,1) - A(2,1)*A(3,4)*A(4,3) - A(2,3)*A(3,1)*A(4,4)) &
          + A(1,3)*(A(2,1)*A(3,2)*A(4,4) + A(2,2)*A(3,4)*A(4,1) + A(2,4)*A(3,1)*A(4,2) &
              - A(2,4)*A(3,2)*A(4,1) - A(2,1)*A(3,4)*A(4,2) - A(2,2)*A(3,1)*A(4,4)) &
          - A(1,4)*(A(2,1)*A(3,2)*A(4,3) + A(2,2)*A(3,3)*A(4,1) + A(2,3)*A(3,1)*A(4,2) &
              - A(2,3)*A(3,2)*A(4,1) - A(2,1)*A(3,3)*A(4,2) - A(2,2)*A(3,1)*A(4,3))

    IF (DET == 0.0) THEN
        PRINT *, "The matrix is singular and cannot be inverted."
        RETURN
    END IF

    ! Compute the cofactors of the 4x4 matrix
    COFACTOR(1,1) = +(A(2,2)*A(3,3)*A(4,4) + A(2,3)*A(3,4)*A(4,2) + A(2,4)*A(3,2)*A(4,3) &
                   - A(2,4)*A(3,3)*A(4,2) - A(2,2)*A(3,4)*A(4,3) - A(2,3)*A(3,2)*A(4,4))
    COFACTOR(1,2) = -(A(2,1)*A(3,3)*A(4,4) + A(2,3)*A(3,4)*A(4,1) + A(2,4)*A(3,1)*A(4,3) &
                   - A(2,4)*A(3,3)*A(4,1) - A(2,1)*A(3,4)*A(4,3) - A(2,3)*A(3,1)*A(4,4))
    COFACTOR(1,3) = +(A(2,1)*A(3,2)*A(4,4) + A(2,2)*A(3,4)*A(4,1) + A(2,4)*A(3,1)*A(4,2) &
                   - A(2,4)*A(3,2)*A(4,1) - A(2,1)*A(3,4)*A(4,2) - A(2,2)*A(3,1)*A(4,4))
    COFACTOR(1,4) = -(A(2,1)*A(3,2)*A(4,3) + A(2,2)*A(3,3)*A(4,1) + A(2,3)*A(3,1)*A(4,2) &
                   - A(2,3)*A(3,2)*A(4,1) - A(2,1)*A(3,3)*A(4,2) - A(2,2)*A(3,1)*A(4,3))

    COFACTOR(2,1) = -(A(1,2)*A(3,3)*A(4,4) + A(1,3)*A(3,4)*A(4,2) + A(1,4)*A(3,2)*A(4,3) &
                   - A(1,4)*A(3,3)*A(4,2) - A(1,2)*A(3,4)*A(4,3) - A(1,3)*A(3,2)*A(4,4))
    COFACTOR(2,2) = +(A(1,1)*A(3,3)*A(4,4) + A(1,3)*A(3,4)*A(4,1) + A(1,4)*A(3,1)*A(4,3) &
                   - A(1,4)*A(3,3)*A(4,1) - A(1,1)*A(3,4)*A(4,3) - A(1,3)*A(3,1)*A(4,4))
    COFACTOR(2,3) = -(A(1,1)*A(3,2)*A(4,4) + A(1,2)*A(3,4)*A(4,1) + A(1,4)*A(3,1)*A(4,2) &
                   - A(1,4)*A(3,2)*A(4,1) - A(1,1)*A(3,4)*A(4,2) - A(1,2)*A(3,1)*A(4,4))
    COFACTOR(2,4) = +(A(1,1)*A(3,2)*A(4,3) + A(1,2)*A(3,3)*A(4,1) + A(1,3)*A(3,1)*A(4,2) &
                   - A(1,3)*A(3,2)*A(4,1) - A(1,1)*A(3,3)*A(4,2) - A(1,2)*A(3,1)*A(4,3))

    COFACTOR(3,1) = +(A(1,2)*A(2,3)*A(4,4) + A(1,3)*A(2,4)*A(4,2) + A(1,4)*A(2,2)*A(4,3) &
                   - A(1,4)*A(2,3)*A(4,2) - A(1,2)*A(2,4)*A(4,3) - A(1,3)*A(2,2)*A(4,4))
    COFACTOR(3,2) = -(A(1,1)*A(2,3)*A(4,4) + A(1,3)*A(2,4)*A(4,1) + A(1,4)*A(2,1)*A(4,3) &
                   - A(1,4)*A(2,3)*A(4,1) - A(1,1)*A(2,4)*A(4,3) - A(1,3)*A(2,1)*A(4,4))
    COFACTOR(3,3) = +(A(1,1)*A(2,2)*A(4,4) + A(1,2)*A(2,4)*A(4,1) + A(1,4)*A(2,1)*A(4,2) &
                   - A(1,4)*A(2,2)*A(4,1) - A(1,1)*A(2,4)*A(4,2) - A(1,2)*A(2,1)*A(4,4))
    COFACTOR(3,4) = -(A(1,1)*A(2,2)*A(4,3) + A(1,2)*A(2,3)*A(4,1) + A(1,3)*A(2,1)*A(4,2) &
                   - A(1,3)*A(2,2)*A(4,1) - A(1,1)*A(2,3)*A(4,2) - A(1,2)*A(2,1)*A(4,3))

    COFACTOR(4,1) = -(A(1,2)*A(2,3)*A(3,4) + A(1,3)*A(2,4)*A(3,2) + A(1,4)*A(2,2)*A(3,3) &
                   - A(1,4)*A(2,3)*A(3,2) - A(1,2)*A(2,4)*A(3,3) - A(1,3)*A(2,2)*A(3,4))
    COFACTOR(4,2) = +(A(1,1)*A(2,3)*A(3,4) + A(1,3)*A(2,4)*A(3,1) + A(1,4)*A(2,1)*A(3,3) &
                   - A(1,4)*A(2,3)*A(3,1) - A(1,1)*A(2,4)*A(3,3) - A(1,3)*A(2,1)*A(3,4))
    COFACTOR(4,3) = -(A(1,1)*A(2,2)*A(3,4) + A(1,2)*A(2,4)*A(3,1) + A(1,4)*A(2,1)*A(3,2) &
                   - A(1,4)*A(2,2)*A(3,1) - A(1,1)*A(2,4)*A(3,2) - A(1,2)*A(2,1)*A(3,4))
    COFACTOR(4,4) = +(A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) &
                   - A(1,3)*A(2,2)*A(3,1) - A(1,1)*A(2,3)*A(3,2) - A(1,2)*A(2,1)*A(3,3))

    ! Transpose the cofactor matrix and divide by the determinant to get the inverse
    AINV = TRANSPOSE(COFACTOR) / DET

  END SUBROUTINE M44INV


!-------------------------------------------------------------------------------!
END SUBROUTINE LeftEigenvectorsObtainedFromRotationalInvariance
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifndef MULTIFLUID
SUBROUTINE RightEigenvectorsXDirection(Cons,RMat)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: nDims
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: Cons(1:nVar)
REAL,INTENT(INOUT) :: RMat(1:nVar,1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: Prim(1:nVar)
REAL               :: ro, u, v, p, Hs, H1, a, abis

#if defined(RELAXATION)
! PRINT*, "Stopping in RightEigenvectorsXDirection. I do not know what to do in IMEX case"
#else

CALL ConsToPrim(Cons(1:nVar),Prim(1:nVar))

ro = Prim(1)
u  = Prim(2)
v  = Prim(3)
p  = Prim(4)
Hs = Gmm*p/(ro*(Gmm-1.0))+0.50*(u**2+v**2)
H1 = -2.0*Hs+u**2+v**2
a  = SQRT(ABS(Gmm*p/ro))




RMat(1:4,1:4)=RESHAPE([1.0, u-a,     v,            Hs-u*a, &
          &   1.0,   u,     v,  0.50*(u**2+v**2), &
          &   0.0, 0.0,   1.0,                 v, &
          &   1.0, u+a,     v,           Hs+u*a],[4,4])


#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE RightEigenvectorsXDirection
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifndef MULTIFLUID
SUBROUTINE LeftEigenvectorsXDirection(Cons,LMat)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: nDims
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: Cons(1:nVar)
REAL,INTENT(INOUT) :: LMat(1:nVar,1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: Prim(1:nVar)
REAL               :: ro, u, v, p, Hs, H1, a, abis

#if defined(RELAXATION)
! PRINT*, "Stopping in LeftEigenvectorsXDirection. I do not know what to do in IMEX case"
#else

CALL ConsToPrim(Cons(1:nVar),Prim(1:nVar))

ro = Prim(1)
u  = Prim(2)
v  = Prim(3)
p  = Prim(4)
Hs = Gmm*p/(ro*(Gmm-1.0))+0.50*(u**2+v**2)
H1 = -2.0*Hs+u**2+v**2
a  = SQRT(ABS(Gmm*p/ro))



LMat(1:4,1:4)=RESHAPE([0.50*(u/a-1.0-2.0*Hs/H1),2.0+2.0*Hs/H1,  -v,-(u+a)/(2.0*a)-Hs/H1, &
            &        -1.0/(2.0*a)+u/H1,    -2.0*u/H1, 0.0,    1.0/(2.0*a)+u/H1, &
            &                     v/H1,    -2.0*v/H1, 1.0,                v/H1, &
            &                  -1.0/H1,       2.0/H1, 0.0,            -1.0/H1],[4,4])

#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE LeftEigenvectorsXDirection
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifndef MULTIFLUID
SUBROUTINE RightEigenvectorsYDirection(Cons,RMat)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: nDims
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: Cons(1:nVar)
REAL,INTENT(INOUT) :: RMat(1:nVar,1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: Prim(1:nVar)
REAL               :: ro, u, v, p, Hs, H1, a, abis

#if defined(RELAXATION)
! PRINT*, "Stopping in RightEigenvectorsYDirection. I do not know what to do in IMEX case"
#else

CALL ConsToPrim(Cons(1:nVar),Prim(1:nVar))

ro = Prim(1)
u  = Prim(2)
v  = Prim(3)
p  = Prim(4)
Hs = Gmm*p/(ro*(Gmm-1.0))+0.50*(u**2+v**2)
H1 = -2.0*Hs+u**2+v**2
a  = SQRT(ABS(Gmm*p/ro))


RMat(1:4,1:4)=RESHAPE([1.0,   u,   v-a,           Hs-v*a, &
            & 1.0,   u,     v, 0.50*(u**2+v**2), &
            & 0.0, 1.0,   0.0,                u, &
            & 1.0,   u,   v+a,          Hs+v*a],[4,4])

#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE RightEigenvectorsYDirection
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifndef MULTIFLUID
SUBROUTINE LeftEigenvectorsYDirection(Cons,LMat)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: nDims
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: Cons(1:nVar)
REAL,INTENT(INOUT) :: LMat(1:nVar,1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: Prim(1:nVar)
REAL               :: ro, u, v, p, Hs, H1, a, abis

#if defined(RELAXATION)
! PRINT*, "Stopping in LeftEigenvectorsYDirection. I do not know what to do in IMEX case"
#else

CALL ConsToPrim(Cons(1:nVar),Prim(1:nVar))

ro = Prim(1)
u  = Prim(2)
v  = Prim(3)
p  = Prim(4)
Hs = Gmm*p/(ro*(Gmm-1.0))+0.50*(u**2+v**2)
H1 = -2.0*Hs+u**2+v**2
a  = SQRT(ABS(Gmm*p/ro))



LMat(1:4,1:4)=RESHAPE([ 0.50*(v/a-1.0-2.0*Hs/H1),    2.0+2.0*Hs/H1,  -u,   -(v+a)/(2.0*a)-Hs/H1, &
              &                    u/H1,        -2.0*u/H1, 1.0,                   u/H1, &
              &       -1.0/(2.0*a)+v/H1,        -2.0*v/H1, 0.0,       1.0/(2.0*a)+v/H1, &
              &                 -1.0/H1,           2.0/H1, 0.0,               -1.0/H1],[4,4])  


#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE LeftEigenvectorsYDirection
#endif
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE RightEigenvectorsRotationalInvariancePrimitiveSystem(Prim,NormVect,RMat)
!*------------------------------------------------------------------------------!
!*Toro's book
!*Proposition 3.15 (Rotational Invariance)
!*K(U,theta)=T^{-1}(theta)K(T(theta)U)
!*where K is the matrix of the eigenvectors of the F flux (first component) of the Euler equations and T is the rotation matrix
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: nDims
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: Prim(1:nVar)
REAL,INTENT(IN)    :: NormVect(nDims)
REAL,INTENT(INOUT) :: RMat(1:nVar,1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: TMat(1:nVar,1:nVar)
REAL               :: invTMat(1:nVar,1:nVar)
REAL               :: KMat(1:nVar,1:nVar)
REAL               :: TU(1:nVar)
REAL               :: theta

#if defined(RELAXATION)
! PRINT*, "Stopping in RightEigenvectorsRotationalInvariancePrimitiveSystem. I do not know what to do in IMEX case"
#else

theta = ATAN2(NormVect(2), NormVect(1))

TMat(1:4,1:4)=RESHAPE([1.0,          0.0,          0.0,               0.0, &
          &   0.0,   COS(theta),  -SIN(theta),               0.0, &
          &   0.0,   SIN(theta),   COS(theta),               0.0, &
          &   0.0,          0.0,          0.0,           1.0],[4,4])

invTMat(1:4,1:4)=RESHAPE([1.0,          0.0,          0.0,               0.0, &
          &   0.0,   COS(theta),   SIN(theta),               0.0, &
          &   0.0,  -SIN(theta),   COS(theta),               0.0, &
          &   0.0,          0.0,          0.0,           1.0],[4,4])


#ifdef MULTIFLUID
TMat(nVar,nVar)    = 1.0
invTMat(nVar,nVar) = 1.0
#endif

TU(1:nVar)=MATMUL(TMat(1:nVar,1:nVar),Prim(1:nVar))

CALL RightEigenvectorsXDirectionPrimitiveSystem(TU(1:nVar),KMat(1:nVar,1:nVar))

RMat(1:nVar,1:nVar)=MATMUL(invTMat(1:nVar,1:nVar),KMat(1:nVar,1:nVar))

#endif


!-------------------------------------------------------------------------------!
END SUBROUTINE RightEigenvectorsRotationalInvariancePrimitiveSystem
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE LeftEigenvectorsRotationalInvariancePrimitiveSystem(Prim,NormVect,LMat)
!*------------------------------------------------------------------------------!
!*Toro's book
!*Proposition 3.15 (Rotational Invariance)
!*K(U,theta)=T^{-1}(theta)K(T(theta)U)
!*where K is the matrix of the eigenvectors of the F flux (first component) of the Euler equations and T is the rotation matrix
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: nDims
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: Prim(1:nVar)
REAL,INTENT(IN)    :: NormVect(1:nDims)
REAL,INTENT(INOUT) :: LMat(1:nVar,1:nVar)
!-------------------------------------------------------------------------------!
REAL               :: RMat(1:nVar,1:nVar)

#if defined(RELAXATION)
! PRINT*, "Stopping in LeftEigenvectorsRotationalInvariancePrimitiveSystem. I do not know what to do in IMEX case"
#else

CALL RightEigenvectorsRotationalInvariancePrimitiveSystem(Prim,NormVect,RMat)
CALL M44INV(RMat,LMat)

#endif



CONTAINS

  SUBROUTINE M44INV (A, AINV)
    IMPLICIT NONE
    REAL, DIMENSION(4,4), INTENT(IN)  :: A
    REAL, DIMENSION(4,4), INTENT(OUT) :: AINV

    REAL :: DET
    REAL, DIMENSION(4,4) :: COFACTOR

    ! Compute the determinant of the 4x4 matrix
    DET =   A(1,1)*(A(2,2)*A(3,3)*A(4,4) + A(2,3)*A(3,4)*A(4,2) + A(2,4)*A(3,2)*A(4,3) &
              - A(2,4)*A(3,3)*A(4,2) - A(2,2)*A(3,4)*A(4,3) - A(2,3)*A(3,2)*A(4,4)) &
          - A(1,2)*(A(2,1)*A(3,3)*A(4,4) + A(2,3)*A(3,4)*A(4,1) + A(2,4)*A(3,1)*A(4,3) &
              - A(2,4)*A(3,3)*A(4,1) - A(2,1)*A(3,4)*A(4,3) - A(2,3)*A(3,1)*A(4,4)) &
          + A(1,3)*(A(2,1)*A(3,2)*A(4,4) + A(2,2)*A(3,4)*A(4,1) + A(2,4)*A(3,1)*A(4,2) &
              - A(2,4)*A(3,2)*A(4,1) - A(2,1)*A(3,4)*A(4,2) - A(2,2)*A(3,1)*A(4,4)) &
          - A(1,4)*(A(2,1)*A(3,2)*A(4,3) + A(2,2)*A(3,3)*A(4,1) + A(2,3)*A(3,1)*A(4,2) &
              - A(2,3)*A(3,2)*A(4,1) - A(2,1)*A(3,3)*A(4,2) - A(2,2)*A(3,1)*A(4,3))

    IF (DET == 0.0) THEN
        PRINT *, "The matrix is singular and cannot be inverted."
        RETURN
    END IF

    ! Compute the cofactors of the 4x4 matrix
    COFACTOR(1,1) = +(A(2,2)*A(3,3)*A(4,4) + A(2,3)*A(3,4)*A(4,2) + A(2,4)*A(3,2)*A(4,3) &
                   - A(2,4)*A(3,3)*A(4,2) - A(2,2)*A(3,4)*A(4,3) - A(2,3)*A(3,2)*A(4,4))
    COFACTOR(1,2) = -(A(2,1)*A(3,3)*A(4,4) + A(2,3)*A(3,4)*A(4,1) + A(2,4)*A(3,1)*A(4,3) &
                   - A(2,4)*A(3,3)*A(4,1) - A(2,1)*A(3,4)*A(4,3) - A(2,3)*A(3,1)*A(4,4))
    COFACTOR(1,3) = +(A(2,1)*A(3,2)*A(4,4) + A(2,2)*A(3,4)*A(4,1) + A(2,4)*A(3,1)*A(4,2) &
                   - A(2,4)*A(3,2)*A(4,1) - A(2,1)*A(3,4)*A(4,2) - A(2,2)*A(3,1)*A(4,4))
    COFACTOR(1,4) = -(A(2,1)*A(3,2)*A(4,3) + A(2,2)*A(3,3)*A(4,1) + A(2,3)*A(3,1)*A(4,2) &
                   - A(2,3)*A(3,2)*A(4,1) - A(2,1)*A(3,3)*A(4,2) - A(2,2)*A(3,1)*A(4,3))

    COFACTOR(2,1) = -(A(1,2)*A(3,3)*A(4,4) + A(1,3)*A(3,4)*A(4,2) + A(1,4)*A(3,2)*A(4,3) &
                   - A(1,4)*A(3,3)*A(4,2) - A(1,2)*A(3,4)*A(4,3) - A(1,3)*A(3,2)*A(4,4))
    COFACTOR(2,2) = +(A(1,1)*A(3,3)*A(4,4) + A(1,3)*A(3,4)*A(4,1) + A(1,4)*A(3,1)*A(4,3) &
                   - A(1,4)*A(3,3)*A(4,1) - A(1,1)*A(3,4)*A(4,3) - A(1,3)*A(3,1)*A(4,4))
    COFACTOR(2,3) = -(A(1,1)*A(3,2)*A(4,4) + A(1,2)*A(3,4)*A(4,1) + A(1,4)*A(3,1)*A(4,2) &
                   - A(1,4)*A(3,2)*A(4,1) - A(1,1)*A(3,4)*A(4,2) - A(1,2)*A(3,1)*A(4,4))
    COFACTOR(2,4) = +(A(1,1)*A(3,2)*A(4,3) + A(1,2)*A(3,3)*A(4,1) + A(1,3)*A(3,1)*A(4,2) &
                   - A(1,3)*A(3,2)*A(4,1) - A(1,1)*A(3,3)*A(4,2) - A(1,2)*A(3,1)*A(4,3))

    COFACTOR(3,1) = +(A(1,2)*A(2,3)*A(4,4) + A(1,3)*A(2,4)*A(4,2) + A(1,4)*A(2,2)*A(4,3) &
                   - A(1,4)*A(2,3)*A(4,2) - A(1,2)*A(2,4)*A(4,3) - A(1,3)*A(2,2)*A(4,4))
    COFACTOR(3,2) = -(A(1,1)*A(2,3)*A(4,4) + A(1,3)*A(2,4)*A(4,1) + A(1,4)*A(2,1)*A(4,3) &
                   - A(1,4)*A(2,3)*A(4,1) - A(1,1)*A(2,4)*A(4,3) - A(1,3)*A(2,1)*A(4,4))
    COFACTOR(3,3) = +(A(1,1)*A(2,2)*A(4,4) + A(1,2)*A(2,4)*A(4,1) + A(1,4)*A(2,1)*A(4,2) &
                   - A(1,4)*A(2,2)*A(4,1) - A(1,1)*A(2,4)*A(4,2) - A(1,2)*A(2,1)*A(4,4))
    COFACTOR(3,4) = -(A(1,1)*A(2,2)*A(4,3) + A(1,2)*A(2,3)*A(4,1) + A(1,3)*A(2,1)*A(4,2) &
                   - A(1,3)*A(2,2)*A(4,1) - A(1,1)*A(2,3)*A(4,2) - A(1,2)*A(2,1)*A(4,3))

    COFACTOR(4,1) = -(A(1,2)*A(2,3)*A(3,4) + A(1,3)*A(2,4)*A(3,2) + A(1,4)*A(2,2)*A(3,3) &
                   - A(1,4)*A(2,3)*A(3,2) - A(1,2)*A(2,4)*A(3,3) - A(1,3)*A(2,2)*A(3,4))
    COFACTOR(4,2) = +(A(1,1)*A(2,3)*A(3,4) + A(1,3)*A(2,4)*A(3,1) + A(1,4)*A(2,1)*A(3,3) &
                   - A(1,4)*A(2,3)*A(3,1) - A(1,1)*A(2,4)*A(3,3) - A(1,3)*A(2,1)*A(3,4))
    COFACTOR(4,3) = -(A(1,1)*A(2,2)*A(3,4) + A(1,2)*A(2,4)*A(3,1) + A(1,4)*A(2,1)*A(3,2) &
                   - A(1,4)*A(2,2)*A(3,1) - A(1,1)*A(2,4)*A(3,2) - A(1,2)*A(2,1)*A(3,4))
    COFACTOR(4,4) = +(A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) &
                   - A(1,3)*A(2,2)*A(3,1) - A(1,1)*A(2,3)*A(3,2) - A(1,2)*A(2,1)*A(3,3))

    ! Transpose the cofactor matrix and divide by the determinant to get the inverse
    AINV = TRANSPOSE(COFACTOR) / DET

  END SUBROUTINE M44INV


!-------------------------------------------------------------------------------!
END SUBROUTINE LeftEigenvectorsRotationalInvariancePrimitiveSystem
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE RightEigenvectorsXDirectionPrimitiveSystem(Prim,RMat)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: nDims
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: Prim(1:nVar)
REAL,INTENT(INOUT) :: RMat(1:nVar,1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: ro, u, v, p, Hs, H1, a, abis, pinf

#if defined(RELAXATION)
! PRINT*, "Stopping in RightEigenvectorsXDirectionPrimitiveSystem. I do not know what to do in IMEX case"
#else

ro = Prim(1)
#ifdef MOMENTUMINPRIMITIVEVARIABLES
u  = Prim(2)/Prim(1)
v  = Prim(3)/Prim(1)
#else
u  = Prim(2)
v  = Prim(3)
#endif
p  = Prim(4)
#ifdef MULTIFLUID
Gmm=Get_Gamma(Prim(5))
pinf=Get_Pressure_Infinity(Prim(5))
a  = SQRT(ABS(Gmm*(p+pinf)/ro))
#else
a  = SQRT(ABS(Gmm*p/ro))
#endif


RMat(1:4,1:4)=RESHAPE([1.0/a**2, -1.0/(ro*a),     0.0,      1.0, &
          &        1.0,         0.0,     0.0,      0.0, &
          &        0.0,         0.0,     1.0,      0.0, &
          &   1.0/a**2,  1.0/(ro*a),     0.0,      1.0],[4,4])


#ifdef MULTIFLUID
RMat(nVar,nVar)=1.0
#endif

#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE RightEigenvectorsXDirectionPrimitiveSystem
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE LeftEigenvectorsXDirectionPrimitiveSystem(Prim,LMat)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: nDims
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: Prim(1:nVar)
REAL,INTENT(INOUT) :: LMat(1:nVar,1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: ro, u, v, p, Hs, H1, a, abis, pinf

#if defined(RELAXATION)
! PRINT*, "Stopping in LeftEigenvectorsXDirectionPrimitiveSystem. I do not know what to do in IMEX case"
#else

ro = Prim(1)
#ifdef MOMENTUMINPRIMITIVEVARIABLES
u  = Prim(2)/Prim(1)
v  = Prim(3)/Prim(1)
#else
u  = Prim(2)
v  = Prim(3)
#endif
p  = Prim(4)
#ifdef MULTIFLUID
Gmm=Get_Gamma(Prim(5))
pinf=Get_Pressure_Infinity(Prim(5))
a  = SQRT(ABS(Gmm*(p+pinf)/ro))
#else
a  = SQRT(ABS(Gmm*p/ro))
#endif


LMat(1:4,1:4)=RESHAPE([                     0.0,          1.0, 0.0,                 0.0, &
            &                -0.5*ro*a,          0.0, 0.0,            0.5*ro*a, &
            &                      0.0,          0.0, 1.0,                 0.0, &
            &                      0.5,    -1.0/a**2, 0.0,                0.5],[4,4])

#ifdef MULTIFLUID
LMat(nVar,nVar)=1.0
#endif

#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE LeftEigenvectorsXDirectionPrimitiveSystem
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE RightEigenvectorsYDirectionPrimitiveSystem(Prim,RMat)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: nDims
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: Prim(1:nVar)
REAL,INTENT(INOUT) :: RMat(1:nVar,1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: ro, u, v, p, Hs, H1, a, abis, pinf

#if defined(RELAXATION)
! PRINT*, "Stopping in RightEigenvectorsYDirectionPrimitiveSystem. I do not know what to do in IMEX case"
#else

ro = Prim(1)
#ifdef MOMENTUMINPRIMITIVEVARIABLES
u  = Prim(2)/Prim(1)
v  = Prim(3)/Prim(1)
#else
u  = Prim(2)
v  = Prim(3)
#endif
p  = Prim(4)
#ifdef MULTIFLUID
Gmm=Get_Gamma(Prim(5))
pinf=Get_Pressure_Infinity(Prim(5))
a  = SQRT(ABS(Gmm*(p+pinf)/ro))
#else
a  = SQRT(ABS(Gmm*p/ro))
#endif

RMat(1:4,1:4)=RESHAPE([1.0/a**2,         0.0,  -1.0/(ro*a),      1.0, &
          &        1.0,         0.0,          0.0,      0.0, &
          &        0.0,        -1.0,          0.0,      0.0, &
          &   1.0/a**2,         0.0,   1.0/(ro*a),      1.0],[4,4])

#ifdef MULTIFLUID
RMat(nVar,nVar)=1.0
#endif

#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE RightEigenvectorsYDirectionPrimitiveSystem
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE LeftEigenvectorsYDirectionPrimitiveSystem(Prim,LMat)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: nDims
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: Prim(1:nVar)
REAL,INTENT(INOUT) :: LMat(1:nVar,1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: ro, u, v, p, Hs, H1, a, abis, pinf

#if defined(RELAXATION)
! PRINT*, "Stopping in LeftEigenvectorsYDirectionPrimitiveSystem. I do not know what to do in IMEX case"
#else

ro = Prim(1)
#ifdef MOMENTUMINPRIMITIVEVARIABLES
u  = Prim(2)/Prim(1)
v  = Prim(3)/Prim(1)
#else
u  = Prim(2)
v  = Prim(3)
#endif
p  = Prim(4)
#ifdef MULTIFLUID
Gmm=Get_Gamma(Prim(5))
pinf=Get_Pressure_Infinity(Prim(5))
a  = SQRT(ABS(Gmm*(p+pinf)/ro))
#else
a  = SQRT(ABS(Gmm*p/ro))
#endif

LMat(1:4,1:4)=RESHAPE([                     0.0,          1.0, 0.0,                 0.0, &
            &                      0.0,          0.0,-1.0,                 0.0, &
            &                -0.5*ro*a,          0.0, 0.0,            0.5*ro*a, &
            &                      0.5,    -1.0/a**2, 0.0,                0.5],[4,4])

#ifdef MULTIFLUID
LMat(nVar,nVar)=1.0
#endif

#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE LeftEigenvectorsYDirectionPrimitiveSystem
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE TransitionMatrixConsToPrim(Cons,LMat)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DENSITY
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: Cons(1:nVar)
REAL,INTENT(INOUT) :: LMat(1:nVar,1:nVar)
!-------------------------------------------------------------------------------!
REAL             :: ro, rou, rov, E, p, u, v, a, H
REAL             :: RMat(1:nVar,1:nVar), MMat(1:nVar,1:nVar)
REAL             :: Prim(1:nVar)
REAL             :: W(1:nVar)
INTEGER          :: iVar, jVar
REAL             :: ueps, veps, roueps, roveps, ro_safe

ro =Cons(1)
rou=Cons(2)
rov=Cons(3)
E  =Cons(4)
u  =rou/ro
v  =rov/ro
p  =(Gmm-1.0)*(E-0.5*ro*(u**2+v**2))

LMat=0.0

#ifdef RELAXATION
ueps=u*EPS_LM**2
veps=v*EPS_LM**2
roueps=rou*EPS_LM**2
roveps=rov*EPS_LM**2
#else
ueps=u
veps=v
roueps=rou
roveps=rov
#endif

#ifdef TRICKSTOAVOIDDIVISIONSBYZERO
IF ( ro .LT. MIN_DENSITY ) THEN
  ro_safe=MIN_DENSITY
ELSE
  ro_safe=ro
END IF
#else
ro_safe=ro
#endif


#ifdef MOMENTUMINPRIMITIVEVARIABLES
LMat(1:4,1) = (/ 1.0,  0.0,    0.0,    0.0               /)
LMat(1:4,2) = (/ 0.0,  1.0,    0.0,  -0.5*(Gmm-1.0)*ueps /)
LMat(1:4,3) = (/ 0.0,  0.0,    1.0,  -0.5*(Gmm-1.0)*veps /)
LMat(1:4,4) = (/ 0.0,  0.0,    0.0,      Gmm-1.0         /)
#else
LMat(1:4,1) = (/ 1.0,         0.0,          0.0,    0.0               /)
LMat(1:4,2) = (/ 0.0, 1.0/ro_safe,          0.0,  -0.5*(Gmm-1.0)*ueps /)
LMat(1:4,3) = (/ 0.0,          0.0, 1.0/ro_safe,  -0.5*(Gmm-1.0)*veps /)
LMat(1:4,4) = (/ 0.0,          0.0,         0.0,      Gmm-1.0         /)
#endif

#ifdef MULTIFLUID
#ifdef PRIMITIVEFORMULATIONLEVELSET
      LMat(5,5)=1.0/ro_safe !*Primitive is phi in this case
#else
      LMat(5,5)=1.0
#endif
#endif


#if(1==0)
PRINT*, "Check TransitionMatrixConsToPrim"
CALL ConsToPrim(Cons(1:nVar),Prim(1:nVar))
W(1:nVar)=MATMUL(LMat(1:nVar,1:nVar),Cons(1:nVar))
DO iVar=1,nVar
  IF( ABS(Prim(iVar)-W(iVar)) .GT. 1e-14 ) THEN
    PRINT*, "PROBLEM", iVar, Prim(iVar), W(iVar)
    STOP
  END IF
END DO
CALL TransitionMatrixPrimToCons(Cons,RMat)
MMat=MATMUL(RMat,LMat)
DO iVar=1,nVar
  DO jVar=1,nVar
    IF (iVar .EQ. jVar) THEN
      IF( ABS(Mmat(iVar,jVar)-1.0) .GT. 1e-14 ) THEN
        PRINT*, "PROBLEM", iVar, jVar, Mmat(iVar,jVar)
        STOP
      END IF
    ELSE
      IF( ABS(Mmat(iVar,jVar)-0.0) .GT. 1e-14 ) THEN
        PRINT*, "PROBLEM", iVar, jVar, Mmat(iVar,jVar)
        STOP
      END IF
    END IF
  END DO
END DO
#endif


!-------------------------------------------------------------------------------!
END SUBROUTINE TransitionMatrixConsToPrim
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE TransitionMatrixPrimToCons(Cons,RMat)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: Cons(1:nVar)
REAL,INTENT(INOUT) :: RMat(1:nVar,1:nVar)
!-------------------------------------------------------------------------------!
REAL             :: ro, rou, rov, E, p, u, v, a, H
REAL             :: LMat(1:nVar,1:nVar), MMat(1:nVar,1:nVar)
REAL             :: Prim(1:nVar)
REAL             :: W(1:nVar)
INTEGER          :: iVar, jVar
REAL             :: ueps, veps, roueps, roveps

ro =Cons(1)
rou=Cons(2)
rov=Cons(3)
E  =Cons(4)
u  =rou/ro
v  =rov/ro
p  =(Gmm-1.0)*(E-0.5*ro*(u**2+v**2))

#ifdef RELAXATION
ueps=u*EPS_LM**2
veps=v*EPS_LM**2
roueps=rou*EPS_LM**2
roveps=rov*EPS_LM**2
#else
ueps=u
veps=v
roueps=rou
roveps=rov
#endif

RMat=0.0

#ifdef MOMENTUMINPRIMITIVEVARIABLES
RMat(1:4,1) = (/ 1.0, 0.0,   0.0,    0.0              /)
RMat(1:4,2) = (/ 0.0, 1.0,   0.0,    0.5*ueps         /)
RMat(1:4,3) = (/ 0.0, 0.0,   1.0,    0.5*veps         /)
RMat(1:4,4) = (/ 0.0, 0.0,   0.0,    1.0/(Gmm-1.0)    /)
#else
RMat(1:4,1) = (/ 1.0, 0.0,   0.0,    0.0              /)
RMat(1:4,2) = (/ 0.0,  ro,   0.0,    0.5*roueps       /)
RMat(1:4,3) = (/ 0.0, 0.0,    ro,    0.5*roveps       /)
RMat(1:4,4) = (/ 0.0, 0.0,   0.0,    1.0/(Gmm-1.0)    /)
#endif

#ifdef MULTIFLUID
#ifdef PRIMITIVEFORMULATIONLEVELSET
      RMat(5,5)=ro !*Primitive is phi in this case
#else
      RMat(5,5)=1.0
#endif
#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE TransitionMatrixPrimToCons
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE NormalFlux1D(Prim,NormVect,TangVect,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#ifdef SW
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
#endif

!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar,1:nGPs)
REAL,INTENT(IN)  :: NormVect(1:nDims)
REAL,INTENT(IN)  :: TangVect(1:nDims)
REAL,INTENT(OUT) :: Flux(1:nVar,1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: PrimRotated(1:nVar,1:nGPs)
INTEGER          :: iGP
!-------------------------------------------------------------------------------!


DO iGP=1,nGPs
  ! Rotating the vector quantities       !
  PrimRotated(1,iGP) = Prim(1,iGP)
  PrimRotated(2,iGP) = NormVect(1)*Prim(2,iGP) + NormVect(2)*Prim(3,iGP)
  PrimRotated(3,iGP) = TangVect(1)*Prim(2,iGP) + TangVect(2)*Prim(3,iGP)
  PrimRotated(4,iGP) = Prim(4,iGP)
#ifdef MULTIFLUID
  PrimRotated(5,iGP) = Prim(5,iGP)
#endif

  CALL EvaluateFlux1D(PrimRotated(1:nVar,iGP),Flux(1:nVar,iGP))

  ! Rotating back the momentum components
  Flux(2:3,iGP) = NormVect(1:nDims)*Flux(2,iGP) &
                + TangVect(1:nDims)*Flux(3,iGP)

END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE NormalFlux1D
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE NumericalFluxPrimitiveSystem(PrimL,PrimR,NormVect,TangVect,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#ifdef SW
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
#endif
USE MOD_FiniteVolume2D_vars,ONLY: WhichRiemannSolverPrimitiveSystem
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
#ifdef MULTIFLUID
REAL,INTENT(IN)  :: PrimL(1:nVar+1,1:nGPs)
REAL,INTENT(IN)  :: PrimR(1:nVar+1,1:nGPs)
#else
REAL,INTENT(IN)  :: PrimL(1:nVar,1:nGPs)
REAL,INTENT(IN)  :: PrimR(1:nVar,1:nGPs)
#endif
REAL,INTENT(IN)  :: NormVect(1:nDims)
REAL,INTENT(IN)  :: TangVect(1:nDims)
REAL,INTENT(OUT) :: Flux(1:nVar,1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: PrimLL(1:nVar,1:nGPs), PrimRR(1:nVar,1:nGPs)
INTEGER          :: iGP
REAL             :: FluidL, FluidR
!-------------------------------------------------------------------------------!



DO iGP=1,nGPs
  ! Rotating the vector quantities       !
  PrimLL(1,iGP) = PrimL(1,iGP)
  PrimLL(2,iGP) = NormVect(1)*PrimL(2,iGP) + NormVect(2)*PrimL(3,iGP)
  PrimLL(3,iGP) = TangVect(1)*PrimL(2,iGP) + TangVect(2)*PrimL(3,iGP)
  PrimLL(4,iGP) = PrimL(4,iGP)
#ifdef MULTIFLUID
  PrimLL(5,iGP) = PrimL(5,iGP)
  FluidL=PrimL(6,iGP)
#endif

  PrimRR(1,iGP) = PrimR(1,iGP)
  PrimRR(2,iGP) = NormVect(1)*PrimR(2,iGP) + NormVect(2)*PrimR(3,iGP)
  PrimRR(3,iGP) = TangVect(1)*PrimR(2,iGP) + TangVect(2)*PrimR(3,iGP)
  PrimRR(4,iGP) = PrimR(4,iGP)
#ifdef MULTIFLUID
  PrimRR(5,iGP) = PrimR(5,iGP)
  FluidR=PrimR(6,iGP)
#endif

  SELECT CASE(WhichRiemannSolverPrimitiveSystem)
    CASE(0)
      CALL RusanovFluxPrimitiveSystem(&
            PrimLL(1:nVar,iGP),PrimRR(1:nVar,iGP),Flux(1:nVar,iGP))
    CASE(-1)
#ifdef MULTIFLUID
      CALL CentralUpwindFluxPrimitiveSystem(&
            PrimLL(1:nVar,iGP),PrimRR(1:nVar,iGP),Flux(1:nVar,iGP),statusL=FluidL,statusR=FluidR)
#else
      CALL CentralUpwindFluxPrimitiveSystem(&
            PrimLL(1:nVar,iGP),PrimRR(1:nVar,iGP),Flux(1:nVar,iGP))
#endif
    CASE(1)
#ifdef MULTIFLUID
      CALL LowDissipationCentralUpwindFluxPrimitiveSystem(&
            PrimLL(1:nVar,iGP),PrimRR(1:nVar,iGP),Flux(1:nVar,iGP),statusL=FluidL,statusR=FluidR)
#else
      CALL LowDissipationCentralUpwindFluxPrimitiveSystem(&
            PrimLL(1:nVar,iGP),PrimRR(1:nVar,iGP),Flux(1:nVar,iGP))
#endif
    CASE DEFAULT
      PRINT*, "Problem in NumericalFluxPrimitiveSystem in equation.f90"
      PRINT*, "Wrong upwinding choice", WhichRiemannSolverPrimitiveSystem
      STOP
    END SELECT

  ! Rotating back the momentum components
  Flux(2:3,iGP) = NormVect(1:nDims)*Flux(2,iGP) &
                + TangVect(1:nDims)*Flux(3,iGP)
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE NumericalFluxPrimitiveSystem
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE CentralUpwindFluxPrimitiveSystem(PrimL,PrimR,Flux,statusL,statusR)
!-------------------------------------------------------------------------------!
USE MOD_Reconstruction,     ONLY: MINMOD
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: EPS_AVOIDING_DIVISION_BY_ZERO_IN_CU
#ifdef IMEX
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: HyperbolicityLoss
#endif
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)          :: PrimL(1:nVar), PrimR(1:nVar)
REAL,INTENT(OUT)         :: Flux(1:nVar)
REAL,INTENT(IN),OPTIONAL :: statusL, statusR
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: FluxL(1:nVar), FluxR(1:nVar)
REAL             :: ap, am                             !*One-sided local speeds of propagation
REAL             :: LambdaMax
REAL             :: roL, roR, vxL, vxR, pL, pR, cL, cR !*Support structures
REAL             :: slowest_L, slowest_R
REAL             :: fastest_L, fastest_R
INTEGER          :: iVar
#ifdef MULTIFLUID
REAL             :: roPhiL, roPhiR
REAL             :: GmmL,   GmmR
REAL             :: pinfL,  pinfR
REAL             :: pinf
#endif
REAL             :: PrimLL(1:nVar), PrimRR(1:nVar)
!-------------------------------------------------------------------------------!

CALL EvaluateFluxPrimitiveSystem1D(PrimL,FluxL)
CALL EvaluateFluxPrimitiveSystem1D(PrimR,FluxR)

!*Support structures
roL = PrimL(1)
roR = PrimR(1)

#ifdef MOMENTUMINPRIMITIVEVARIABLES
vxL = PrimL(2)/PrimL(1)
vxR = PrimR(2)/PrimR(1)
#else
vxL = PrimL(2)
vxR = PrimR(2)
#endif

pL  = PrimL(4)
pR  = PrimR(4)


PrimLL(1:nVar)=PrimL(1:nVar)
PrimRR(1:nVar)=PrimR(1:nVar)


#ifdef PWLFORGAMMAPINF
PrimLL(nVar)=PrimL(nVar)
PrimRR(nVar)=PrimR(nVar)
#else
#ifdef MULTIFLUID
IF(PRESENT(statusL)) THEN
  roPhiL = statusL
END IF
IF(PRESENT(statusR)) THEN
  roPhiR = statusR
END IF
PrimLL(nVar)=statusL
PrimRR(nVar)=statusR
#endif
#endif

CALL WaveSpeeds1D(PrimLL,slowest=slowest_L,fastest=fastest_L)
CALL WaveSpeeds1D(PrimRR,slowest=slowest_R,fastest=fastest_R)


!*One-sided local speeds of propagation
am = MIN( slowest_L, slowest_R, -EPS_AVOIDING_DIVISION_BY_ZERO_IN_CU )
ap = MAX( fastest_L, fastest_R, +EPS_AVOIDING_DIVISION_BY_ZERO_IN_CU )

Flux = ( ap*FluxL - am*FluxR +ap*am*(PrimR-PrimL) ) / (ap-am)

!-------------------------------------------------------------------------------!
END SUBROUTINE CentralUpwindFluxPrimitiveSystem
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE LowDissipationCentralUpwindFluxPrimitiveSystem(PrimL,PrimR,Flux,statusL,statusR)
!-------------------------------------------------------------------------------!
USE MOD_Reconstruction,     ONLY: MINMOD
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: EPS_AVOIDING_DIVISION_BY_ZERO_IN_CU
#ifdef IMEX
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: HyperbolicityLoss
#endif
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
USE MOD_FiniteVolume2D_vars,ONLY: WhichRiemannSolver
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)          :: PrimL(1:nVar), PrimR(1:nVar)
REAL,INTENT(OUT)         :: Flux(1:nVar)
REAL,INTENT(IN),OPTIONAL :: statusL, statusR
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: FluxL(1:nVar), FluxR(1:nVar)
REAL             :: ConsL(1:nVar), ConsR(1:nVar), Cstar(1:nVar), w1, w2, SL, SR
REAL             :: ap, am                             !*One-sided local speeds of propagation
REAL             :: LambdaMax
REAL             :: roL, roR, vxL, vxR, pL, pR, cL, cR !*Support structures
REAL             :: slowest_L, slowest_R
REAL             :: fastest_L, fastest_R
REAL             :: dV(1:nVar)                         !*Built-in antidiffusion term
REAL             :: Vstar(1:nVar)                      !*Intermediate value of V
INTEGER          :: iVar
#ifdef MULTIFLUID
REAL             :: roPhiL, roPhiR
REAL             :: GmmL,   GmmR
REAL             :: pinfL,  pinfR
REAL             :: pinf
#endif
REAL             :: PrimLL(1:nVar), PrimRR(1:nVar)
REAL             :: Q(1:nVar)                         !*New antidiffusion term
REAL             :: alphastar, qro, uxstar,uystar
!-------------------------------------------------------------------------------!

CALL EvaluateFluxPrimitiveSystem1D(PrimL,FluxL)
CALL EvaluateFluxPrimitiveSystem1D(PrimR,FluxR)

!*Support structures
roL = PrimL(1)
roR = PrimR(1)

#ifdef MOMENTUMINPRIMITIVEVARIABLES
vxL = PrimL(2)/PrimL(1)
vxR = PrimR(2)/PrimR(1)
#else
vxL = PrimL(2)
vxR = PrimR(2)
#endif

pL  = PrimL(4)
pR  = PrimR(4)


PrimLL(1:nVar)=PrimL(1:nVar)
PrimRR(1:nVar)=PrimR(1:nVar)


#ifdef PWLFORGAMMAPINF
PrimLL(nVar)=PrimL(nVar)
PrimRR(nVar)=PrimR(nVar)
#else
#ifdef MULTIFLUID
IF(PRESENT(statusL)) THEN
  roPhiL = statusL
END IF
IF(PRESENT(statusR)) THEN
  roPhiR = statusR
END IF
PrimLL(nVar)=statusL
PrimRR(nVar)=statusR
#endif
#endif

CALL WaveSpeeds1D(PrimLL,slowest=slowest_L,fastest=fastest_L)
CALL WaveSpeeds1D(PrimRR,slowest=slowest_R,fastest=fastest_R)


!*One-sided local speeds of propagation
SL = MIN( slowest_L, slowest_R, -EPS_AVOIDING_DIVISION_BY_ZERO_IN_CU )
SR = MAX( fastest_L, fastest_R, +EPS_AVOIDING_DIVISION_BY_ZERO_IN_CU )
am=SL
ap=SR


#ifdef NEWANTIDIFFUSIONTERM
PRINT*, "DO NOT USE"
STOP
Q           = 0.0
Vstar       = ( ap*PrimR - am*PrimL - ( FluxR - FluxL ) ) / ( ap-am )
Vstar(4)    = 0.0
Vstar(nVar) = 0.0 !*Only active for MULTIFLUID
#ifdef MOMENTUMINPRIMITIVEVARIABLES
uxstar       = Vstar(2)/Vstar(1)
uystar       = Vstar(3)/Vstar(1)
#else
uxstar       = Vstar(2)
uystar       = Vstar(3)
#endif
IF (uxstar .LT. 0.0) THEN
  alphastar=ap/(ap-uxstar)
ELSE
  alphastar=am/(am-uxstar)
END IF
qro=MINMOD( (uxstar-am)*(Vstar(1)-PrimL(1)), (ap-uxstar)*(PrimR(1)-Vstar(1)) )
Q(1)=alphastar*qro*1.0
#ifdef MOMENTUMINPRIMITIVEVARIABLES
Q(2)=alphastar*qro*uxstar
Q(3)=alphastar*qro*uystar
#else
Q(2)=alphastar*uxstar
Q(3)=alphastar*uystar
#endif

Flux = ( ap*FluxL - am*FluxR +ap*am*(PrimR-PrimL) ) / (ap-am) + Q
#else
Vstar = ( ap*PrimR - am*PrimL - ( FluxR - FluxL ) ) / ( ap-am )
dV = 0.0
Do iVar=1,nVar
  dV(iVar)=MINMOD( Vstar(iVar)-PrimL(iVar), PrimR(iVar)-Vstar(iVar) )
END DO
Flux = ( ap*FluxL - am*FluxR +ap*am*(PrimR-PrimL-dV) ) / (ap-am)
#endif

#ifndef OLDEVOLVERHOWITHSAMESCHEME
#ifdef EVOLVERHOWITHSAMESCHEME
IF (WhichRiemannSolver .EQ. 1) THEN
  !*Do nothing
ELSEIF (WhichRiemannSolver .EQ. 2) THEN

  !*Compute new antidiffusion term for first component to make sure the update is done in the same way
  CALL PrimToCons(PrimL(1:nVar),ConsL(1:nVar))
  CALL PrimToCons(PrimR(1:nVar),ConsR(1:nVar))

  !*Built-in antidiffusion term
  CALL EvaluateFlux1D(PrimL,FluxL)
  CALL EvaluateFlux1D(PrimR,FluxR)

  !*Built-in antidiffusion term
  !*=============================
  !*www=axp(j)-axm(j)
  !*rhost=(axp(j)*uW(1,j1)-axm(j)*uE(1,j)-fW(1,j1)+fE(1,j))/www
  !*rhoust=(axp(j)*uW(2,j1)-axm(j)*uE(2,j)-fW(2,j1)+fE(2,j))/www
  !*=============================
  Cstar = ( SR*ConsR - SL*ConsL - ( FluxR - FluxL ) ) / ( SR - SL )

  !*=============================
  !*ust=rhoust/rhost
  !*=============================
  uxstar = Cstar(2)/Cstar(1)
  uystar = Cstar(3)/Cstar(1)

  !*=============================
  !*if(ust<0.d0) then
  !*  alphast=axp(j)/(axp(j)-ust)
  !*else
  !*  alphast=axm(j)/(axm(j)-ust)
  !*endif
  !*=============================
  IF (uxstar .LT. 0.0) THEN
    alphastar=ap/(ap-uxstar)
  ELSE
    alphastar=am/(am-uxstar)
  END IF

  !*=============================
  !*w1=(ust-axm(j))*(rhost-vE(1,j))
  !*w2=(axp(j)-ust)*(vW(1,j1)-rhost)
  !*=============================
  w1  = (uxstar-am)*(Cstar(1)-PrimL(1))
  w2  = (ap-uxstar)*(PrimR(1)-Cstar(1)) 

  !*=============================
  !*q(1)=alphast*dmm(w1,w2)
  !*=============================
  qro = MINMOD(w1,w2)
  Q(1)=alphastar*qro

  !*=============================
  !*q(2)=q(1)*ust
  !*=============================
  Q(2)=Q(1)*uxstar
  Q(3)=Q(1)*uystar

  !*=============================
  !*q(3)=0.5d0*q(2)*ust
  !*=============================
  Q(4)=0.5*(q(2)*uxstar+q(3)*uystar)


  !*=============================
  !*do m=1,3
  !* Hx(m,j)=(axp(j)*fE(m,j)-axm(j)*fW(m,j1)+&
  !*          axp(j)*axm(j)*(uW(m,j1)-uE(m,j)))/www+q(m)
  !*enddo
  !*=============================
  Flux(1) = ( ap*FluxL(1) - am*FluxR(1) +ap*am*(ConsR(1)-ConsL(1)) ) / (ap-am) + Q(1) !*NB: ONLY FIRST COMPONENT

ELSE
  PRINT*, "You should not be here in LowDissipationCentralUpwindFluxPrimitiveSystem"
  PRINT*, "The flag EVOLVERHOWITHSAMESCHEME requires a LowDissipationCentralUpwind also for the conserved system"
  STOP
END IF
#endif
#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE LowDissipationCentralUpwindFluxPrimitiveSystem
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE RusanovFluxPrimitiveSystem(PrimL,PrimR,Flux)
!-------------------------------------------------------------------------------!
USE MOD_Reconstruction,     ONLY: MINMOD
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: EPS_AVOIDING_DIVISION_BY_ZERO_IN_CU
#ifdef IMEX
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: HyperbolicityLoss
#endif
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: PrimL(1:nVar), PrimR(1:nVar)
REAL,INTENT(OUT) :: Flux(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: FluxL(1:nVar), FluxR(1:nVar)
REAL             :: ap, am                             !*One-sided local speeds of propagation
REAL             :: LambdaMax
REAL             :: roL, roR, vxL, vxR, pL, pR, cL, cR !*Support structures
REAL             :: slowest_L, slowest_R
REAL             :: fastest_L, fastest_R
REAL             :: dV(1:nVar)                         !*Built-in antidiffusion term
REAL             :: Vstar(1:nVar)                      !*Intermediate value of V
INTEGER          :: iVar
#ifdef MULTIFLUID
REAL             :: roPhiL, roPhiR
REAL             :: GmmL,   GmmR
REAL             :: pinfL,  pinfR
REAL             :: pinf
#endif
!-------------------------------------------------------------------------------!

CALL EvaluateFluxPrimitiveSystem1D(PrimL,FluxL)
CALL EvaluateFluxPrimitiveSystem1D(PrimR,FluxR)

!*Support structures
roL = PrimL(1)
roR = PrimR(1)

#ifdef MOMENTUMINPRIMITIVEVARIABLES
vxL = PrimL(2)/PrimL(1)
vxR = PrimR(2)/PrimR(1)
#else
vxL = PrimL(2)
vxR = PrimR(2)
#endif

pL  = PrimL(4)
pR  = PrimR(4)

#ifdef MULTIFLUID
roPhiL = PrimL(5)
roPhiR = PrimR(5)
#endif

CALL WaveSpeeds1D(PrimL,slowest=slowest_L,fastest=fastest_L)
CALL WaveSpeeds1D(PrimR,slowest=slowest_R,fastest=fastest_R)


!*One-sided local speeds of propagation
am = MIN( slowest_L, slowest_R, -EPS_AVOIDING_DIVISION_BY_ZERO_IN_CU )
ap = MAX( fastest_L, fastest_R, +EPS_AVOIDING_DIVISION_BY_ZERO_IN_CU )

LambdaMax=MAX(ABS(am),ABS(ap))
Flux = 0.5*((FluxL + FluxR) - LambdaMax*(PrimR - PrimL))

!-------------------------------------------------------------------------------!
END SUBROUTINE RusanovFluxPrimitiveSystem
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE EvaluateFluxPrimitiveSystem1D(Prim,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DENSITY, MIN_SPEED, MIN_PRESSURE
#ifdef SW
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
#endif
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: Flux(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: ro, vx, vy, p!*, Energy
#ifdef MULTIFLUID
REAL             :: roPhi
REAL             :: pinf
#endif
!-------------------------------------------------------------------------------!

ro = Prim(1)
#ifdef MOMENTUMINPRIMITIVEVARIABLES
vx = Prim(2)/Prim(1)
vy = Prim(3)/Prim(1)
#else
vx = Prim(2)
vy = Prim(3)
#endif
p  = Prim(4)
#ifdef SW
p  = Kappa*ro**Gmm
#endif
#ifdef MULTIFLUID
roPhi=Prim(5)
#endif

!*Not needbut just in case
#ifdef MULTIFLUID
Gmm=Get_Gamma(roPhi)
pinf=Get_Pressure_Infinity(roPhi)
#endif




#ifdef PATANKAR

#else
#ifdef POSITIVITYCHECKS
IF (ro .LT. MIN_DENSITY) THEN
 ro = MIN_DENSITY
 vx=0.
 vy=0.
END IF
#endif
#endif


#ifdef POSITIVITYCHECKS
IF (p .LT. MIN_PRESSURE) THEN
 p = MIN_PRESSURE
 vx=0.
 vy=0.
END IF
#endif

!*Energy = p/(Gmm-1.0)+0.5*ro*(vx**2+vy**2)
!*Energy = p/(Gmm-1.0)+EPS_LM**2*0.5*ro*(vx**2+vy**2)

Flux(1) = ro*vx
#ifdef MOMENTUMINPRIMITIVEVARIABLES
!*START MOMENTUMINPRIMITIVEVARIABLES

#ifdef IMEXMOMENTUM
  !*START IMEXMOMENTUM
  Flux(2) = ro*vx**2

#else
  !*ELSE IMEXMOMENTUM

#ifdef RELAXATION
    !*START RELAXATION
    Flux(2) = ro*vx**2+p/EPS_LM**2
#else
    !*ELSE RELAXATION
    Flux(2) = ro*vx**2+p
#endif
    !*END RELAXATION

#endif
  !*END IMEXMOMENTUM

  Flux(3) = ro*vx*vy


#else
!*ELSE MOMENTUMINPRIMITIVEVARIABLES

Flux(2) = 0.5*vx**2
Flux(3) = 0.0

#endif
!*END MOMENTUMINPRIMITIVEVARIABLES

!*=======================
!*SPLITTING
!*=======================
#if defined(ALTERNATIVEFORMULATIONPRESSURE)
!*PRINT*, "This is a consistent formulation."
!*PRINT*, "It works but it is not suitable for the splitting."
!*STOP
Flux(4) = p*vx
#ifdef IMEX
PRINT*, "In particular, for IMEX one really can't do the splitting in this case."
STOP
#endif
#else
Flux(4) = 0.0 
#endif

#ifdef MULTIFLUID
#ifndef PRIMITIVEFORMULATIONLEVELSET
Flux(5) = roPhi*vx
#else
Flux(5) = 0.0 !*Non-conserved product in this case
#endif
#endif


!-------------------------------------------------------------------------------!
END SUBROUTINE EvaluateFluxPrimitiveSystem1D
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE PathConservativeSurfaceContribution(PrimL,PrimR,NormVect,TangVect,B_Sur,am,ap)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#ifdef SW
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
#ifdef MULTIFLUID
REAL,INTENT(IN)  :: PrimL(1:nVar+1,1:nGPs)
REAL,INTENT(IN)  :: PrimR(1:nVar+1,1:nGPs)
#else
REAL,INTENT(IN)  :: PrimL(1:nVar,1:nGPs)
REAL,INTENT(IN)  :: PrimR(1:nVar,1:nGPs)
#endif
REAL,INTENT(IN)  :: NormVect(1:nDims)
REAL,INTENT(IN)  :: TangVect(1:nDims)
REAL,INTENT(OUT) :: B_Sur(1:nVar,1:nGPs)
REAL,INTENT(OUT) :: am(1:nGPs)
REAL,INTENT(OUT) :: ap(1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: PrimLL(1:nVar,1:nGPs), PrimRR(1:nVar,1:nGPs)
INTEGER          :: iGP
REAL             :: FluidL, FluidR
!-------------------------------------------------------------------------------!



DO iGP=1,nGPs
  ! Rotating the vector quantities       !
  PrimLL(1,iGP) = PrimL(1,iGP)
  PrimLL(2,iGP) = NormVect(1)*PrimL(2,iGP) + NormVect(2)*PrimL(3,iGP)
  PrimLL(3,iGP) = TangVect(1)*PrimL(2,iGP) + TangVect(2)*PrimL(3,iGP)
  PrimLL(4,iGP) = PrimL(4,iGP)
#ifdef MULTIFLUID
  PrimLL(5,iGP) = PrimL(5,iGP)
  FluidL=PrimL(6,iGP)
#endif

  PrimRR(1,iGP) = PrimR(1,iGP)
  PrimRR(2,iGP) = NormVect(1)*PrimR(2,iGP) + NormVect(2)*PrimR(3,iGP)
  PrimRR(3,iGP) = TangVect(1)*PrimR(2,iGP) + TangVect(2)*PrimR(3,iGP)
  PrimRR(4,iGP) = PrimR(4,iGP)
#ifdef MULTIFLUID
  PrimRR(5,iGP) = PrimR(5,iGP)
  FluidR=PrimR(6,iGP)
#endif

#ifdef MULTIFLUID
  CALL PathConservativeLinearPath(PrimLL(1:nVar,iGP),PrimRR(1:nVar,iGP),B_Sur(1:nVar,iGP),am(iGP),ap(iGP),statusL=FluidL,statusR=fluidR)
#else
  CALL PathConservativeLinearPath(PrimLL(1:nVar,iGP),PrimRR(1:nVar,iGP),B_Sur(1:nVar,iGP),am(iGP),ap(iGP))
#endif

  ! Rotating back the momentum components
  B_Sur(2:3,iGP) = NormVect(1:nDims)*B_Sur(2,iGP) &
                 + TangVect(1:nDims)*B_Sur(3,iGP)
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE PathConservativeSurfaceContribution
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE PathConservativeLinearPath(PrimL,PrimR,B_Sur,am,ap,statusL,statusR)
!-------------------------------------------------------------------------------!
USE MOD_Reconstruction,     ONLY: MINMOD
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: EPS_AVOIDING_DIVISION_BY_ZERO_IN_CU
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
#ifdef IMEX
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: min_p
USE MOD_FiniteVolume2D_vars,ONLY: HyperbolicityLoss
#endif
USE MOD_FiniteVolume2D_vars,ONLY: timescheme
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)            :: PrimL(1:nVar), PrimR(1:nVar)
REAL,INTENT(OUT)           :: B_Sur(1:nVar)
REAL,INTENT(OUT)           :: am, ap  !*One-sided local speeds of propagation
REAL,INTENT(IN), OPTIONAL  :: statusL, statusR
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: roL, roR, vxL, vxR, vyL, vyR, pL, pR, cL, cR      !*Support structures
REAL             :: slowest_L, slowest_R
REAL             :: fastest_L, fastest_R
INTEGER          :: iVar
#ifdef MULTIFLUID
REAL             :: roPhiL, roPhiR, GmmL, GmmR
REAL             :: pinfL, pinfR
REAL             :: pinf
#endif
REAL             :: PrimLL(1:nVar), PrimRR(1:nVar)
!-------------------------------------------------------------------------------!

B_Sur=0.0

!*Support structures
roL = PrimL(1)
roR = PrimR(1)

#ifdef MOMENTUMINPRIMITIVEVARIABLES
vxL = PrimL(2)/PrimL(1)
vxR = PrimR(2)/PrimR(1)

vyL = PrimL(3)/PrimL(1)
vyR = PrimR(3)/PrimR(1)
#else
vxL = PrimL(2)
vxR = PrimR(2)

vyL = PrimL(3)
vyR = PrimR(3)
#endif

pL  = PrimL(4)
pR  = PrimR(4)

PrimLL(1:nVar)=PrimL(1:nVar)
PrimRR(1:nVar)=PrimR(1:nVar)

#ifdef PWLFORGAMMAPINF
roPhiL = PrimL(nVar) 
roPhiR = PrimR(nVar) 
#else
#ifdef MULTIFLUID
IF(PRESENT(statusL)) THEN
  roPhiL = statusL
END IF
IF(PRESENT(statusR)) THEN
  roPhiR = statusR
END IF
PrimLL(nVar)=statusL
PrimRR(nVar)=statusR
#endif
#endif


#ifdef MULTIFLUID
Gmm =Get_Gamma(roPhiL)
pinf=Get_Pressure_Infinity(roPhiL)
GmmL=Gmm
pinfL=pinf
#endif


#ifdef MULTIFLUID
Gmm=Get_Gamma(roPhiR)
pinf=Get_Pressure_Infinity(roPhiR)
GmmR=Gmm
pinfR=pinf
#endif

CALL WaveSpeeds1D(PrimLL,slowest=slowest_L,fastest=fastest_L)
CALL WaveSpeeds1D(PrimRR,slowest=slowest_R,fastest=fastest_R)


!*One-sided local speeds of propagation
am = MIN( slowest_L, slowest_R, -EPS_AVOIDING_DIVISION_BY_ZERO_IN_CU )
ap = MAX( fastest_L, fastest_R, +EPS_AVOIDING_DIVISION_BY_ZERO_IN_CU )


B_Sur(1) =  0.0


#ifndef MOMENTUMINPRIMITIVEVARIABLES
!*START not MOMENTUMINPRIMITIVEVARIABLES

#ifdef IMEX
  !*START IMEX

#ifdef IMEXL2FULLYUPWINDED
    !*START IMEXL2FULLYUPWINDED
    B_Sur(2) = -0.5 *             ( 1.0/roL + 1.0/roR ) * ( pR  - pL  )
#else
    !*ELSE IMEXL2FULLYUPWINDED
    B_Sur(2) = -0.5 *             ( (max_ro-roL)/(roL*max_ro) + (max_ro-roR)/(roR*max_ro) ) * ( pR  - pL  )
#endif
    !*END IMEXL2FULLYUPWINDED

#else
  !*ELSE IMEX
  B_Sur(2) = -0.5 *             ( 1.0/roL + 1.0/roR ) * ( pR  - pL  )
#endif
  !*END IMEX

#ifdef RELAXATION
  !*START RELAXATION
  B_Sur(2) = B_Sur(2)/EPS_LM**2 !*In case I divide by EPS_LM**2
#endif
  !*END RELAXATION

B_Sur(3) = -0.5 *             ( vxL + vxR         ) * ( vyR - vyL )

#endif
!*END not MOMENTUMINPRIMITIVEVARIABLES


!*=======================
!*SPLITTING
!*=======================
#if defined(ALTERNATIVEFORMULATIONPRESSURE)
!*START ALTERNATIVEFORMULATIONPRESSURE

! PRINT*, "This is a consistent formulation."
! PRINT*, "It works but it is not suitable for the splitting."
! STOP
#ifdef MULTIFLUID
  !*IF MULTIFLUID
  B_Sur(4) =  - 0.5 * ( (GmmL-1.0)*(pL) + (GmmR-1.0)*(pR) + GmmL*pinfL + GmmR*pinfR ) * ( vxR  - vxL  )
#else
  !*ELSE MULTIFLUID
  B_Sur(4) =  - 0.5 * (Gmm-1.0)*( pL + pR         ) * ( vxR  - vxL  )
#endif
  !*END MULTIFLUID

#ifdef IMEX
  !*IF IMEX
  PRINT*, "In particular, for IMEX one really can't do the splitting in this case."
  STOP
#endif
  !*END IMEX




#elif defined(PRESSUREFORMULATIONFORIMEXMOMENTUM)
!*ELSE ALTERNATIVEFORMULATIONPRESSURE

#ifdef IMEXMOMENTUM
    !*START IMEXMOMENTUM
    B_Sur(4) =  + 0.5 * Gmm* ( pL/roL*vxL + pR/roR*vxR ) * ( roR  - roL  ) &
            &  - 0.5 *      (        vxL + vxR        ) * (  pR  - pL   )
#else
    !*ELSE IMEXMOMENTUM
    B_Sur(4) =  + 0.5 * Gmm* ( pL/roL*vxL + pR/roR*vxR ) * ( roR  - roL  ) &
            &  - 0.5 * Gmm* (     pL/roL + pR/roR     ) * ( roR*vxR  - roL*vxL  ) &
            &  - 0.5 *      (        vxL + vxR        ) * (  pR  - pL   )
#endif
    !*END IMEXMOMENTUM

#else
!*ELSE ALTERNATIVEFORMULATIONPRESSURE

#ifdef IMEX
  !*START IMEX

#ifdef IMEXL2FULLYUPWINDED
    !*START IMEXL2FULLYUPWINDED
    B_Sur(4)  =  -0.5 * Gmm * ( pL  + pR  ) * ( vxR - vxL ) & !*This must be splitted
              &  -0.5       * ( vxL + vxR ) * ( pR  - pL  )  
#else
    !*ELSE IMEXL2FULLYUPWINDED
    B_Sur(4)  =  -0.5 * Gmm * ( (pL - min_p)  + (pR - min_p)  ) * ( vxR - vxL ) & !*This must be splitted 
              &  -0.5       * ( vxL + vxR ) * ( pR  - pL  )                       !*NB: This other one not
#endif
    !*END IMEXL2FULLYUPWINDED

#else
  !*ELSE IMEX

#ifdef MULTIFLUID
    !*IF MULTIFLUID
    B_Sur(4)  =  -0.5 * ( GmmL * (pL+pinfL)  + GmmR * (pR+pinfR)  ) * ( vxR - vxL ) & !*This must be splitted
              &  -0.5       * ( vxL + vxR ) * ( pR  - pL  )  

#ifndef PRIMITIVEFORMULATIONLEVELSET
      !*IF not PRIMITIVEFORMULATIONLEVELSET
      B_Sur(nVar)  = 0.0
#else
      !*ELSE not PRIMITIVEFORMULATIONLEVELSET
      B_Sur(nVar)  = - 0.5 * ( vxL + vxR ) * ( roPhiR  - roPhiL  ) !*RMK: roPhi is Phi in this case
#endif
      !*END not PRIMITIVEFORMULATIONLEVELSET

#else
    !*ELSE MULTIFLUID
    B_Sur(4)  =  -0.5 * Gmm * ( pL  + pR  ) * ( vxR - vxL ) & !*This must be splitted
              &  -0.5       * ( vxL + vxR ) * ( pR  - pL  )  
#endif
    !*END MULTIFLUID


#endif
  !*END IMEX

#endif
!*END ALTERNATIVEFORMULATIONPRESSURE

!-------------------------------------------------------------------------------!
END SUBROUTINE PathConservativeLinearPath
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#ifdef IMEX
SUBROUTINE Compute_min_p_max_ro()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
#ifdef ACTIVEFLUX
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
#endif
#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
#endif
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: min_p, max_ro
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
USE MOD_FiniteVolume2D_vars,ONLY: SAFETY_K_TO_AVOID_HYPERBOLICITY_LOSS
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj
REAL    :: Prim(1:nVar)
!-------------------------------------------------------------------------------!

min_p  = HUGE(1.0)
max_ro = 0.0

#if(1==0)
#ifndef PRIMITIVEONLY
!*=======================
!*ALERT: CONSERVED VARIABLES EXCLUDED
!*=======================
!*U
DO jj=1,nElemsY
  DO ii=1,nElemsX
    CALL ConsToPrim(U(1:nVar,ii,jj),Prim(1:nVar))
    max_ro = MAX(max_ro,U(1,ii,jj))
    min_p  = MIN(min_p,Prim(4))
  END DO
END DO
#endif
#endif

#ifdef ACTIVEFLUX
!*W_X
DO jj=1,nElemsY
  DO ii=1,nElemsX+1 !*NB: +1 in X direction
    max_ro = MAX(max_ro,W_X(1,ii,jj))
    min_p  = MIN(min_p,W_X(4,ii,jj))
  END DO
END DO

!*W_Y
DO jj=1,nElemsY+1 !*NB: +1 in Y direction
  DO ii=1,nElemsX
    max_ro = MAX(max_ro,W_Y(1,ii,jj))
    min_p  = MIN(min_p,W_Y(4,ii,jj))
  END DO
END DO
#endif

#ifdef CENTEREDPRIMITIVE
DO jj=1,nElemsY
  DO ii=1,nElemsX
    max_ro = MAX(max_ro,WC(1,ii,jj))
    min_p  = MIN(min_p,WC(4,ii,jj))
  END DO
END DO
#endif



#ifdef NOMODIFICATIONWAVESPEEDIMEX

!*DO NOT TOUCH max_ro AND min_p

#elif HYPERBOLICITYTRICKMODIFIED

#ifdef RELAXATION
max_ro = max_ro*(1.0+K_function(EPS_LM))
min_p  = min_p *(1.0-K_function(EPS_LM))
#else
max_ro = max_ro*(1.0+K_function(1.0))
min_p  = min_p *(1.0-K_function(1.0))
#endif

#else

#ifdef RELAXATION
max_ro = max_ro+EPS_LM*SAFETY_K_TO_AVOID_HYPERBOLICITY_LOSS
min_p  = min_p -EPS_LM*SAFETY_K_TO_AVOID_HYPERBOLICITY_LOSS
#else
max_ro = max_ro+1.0*SAFETY_K_TO_AVOID_HYPERBOLICITY_LOSS
min_p  = min_p -1.0*SAFETY_K_TO_AVOID_HYPERBOLICITY_LOSS
#endif

#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE Compute_min_p_max_ro
#endif
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#ifdef IMEX
FUNCTION K_function(x)
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: x
REAL            :: K_function
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!

K_function = x**4

!-------------------------------------------------------------------------------!
END FUNCTION K_function
#endif
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE Impose_BC_on_Wt()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsType
!*For staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X

!*For staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
INTEGER            :: idx_vx, idx_vy
REAL               :: x0, xc, xt
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

idx_vx = 2
idx_vy = 3

!------------------------------!
! Left Boundary Conditions     !
!------------------------------!
!*Staggering in X direction
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=0,nGhosts 
        Wt_X(1:nVar,-nGhosts+ii,jj) = Wt_X(1:nVar,nElemsX-nGhosts+ii,jj) !*NB: NO CHANGE IN INDICES
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        Wt_X(1:nVar,-nGhosts+ii,jj) = Wt_X(1:nVar,nGhosts-ii+1+1,jj) !*NB: +1 at the RHS
      END DO
    END DO
  CASE(3) ! Inflow
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        Wt_X(1:nVar,-nGhosts+ii,jj) = 0.0
      END DO
    END DO
  CASE(4) ! Outflow
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        Wt_X(1:nVar,-nGhosts+ii,jj) = 0.0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        Wt_X(1:nVar,-nGhosts+ii,jj) = Wt_X(1:nVar,nGhosts-ii+1+1,jj) !*NB:+1 at the RHS
        Wt_X(idx_vx,-nGhosts+ii,jj) =-Wt_X(idx_vx,nGhosts-ii+1+1,jj) !*NB:+1 at the RHS
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition L not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!*Staggering in Y direction
!*NB: +1 in loop in Y direction
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=0,nGhosts
        Wt_Y(1:nVar,-nGhosts+ii,jj) = Wt_Y(1:nVar,nElemsX-nGhosts+ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=0,nGhosts
        Wt_Y(1:nVar,-nGhosts+ii,jj) = Wt_Y(1:nVar,nGhosts-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=0,nGhosts
        Wt_Y(1:nVar,-nGhosts+ii,jj) = 0.0
      END DO
    END DO
  CASE(4) ! Outflow
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=0,nGhosts
        Wt_Y(1:nVar,-nGhosts+ii,jj) = 0.0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=0,nGhosts
        Wt_Y(1:nVar,-nGhosts+ii,jj) = Wt_Y(1:nVar,nGhosts-ii+1,jj)
        Wt_Y(idx_vx,-nGhosts+ii,jj) =-Wt_Y(idx_vx,nGhosts-ii+1,jj)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition L not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!------------------------------!
! Right Boundary Conditions    !
!------------------------------!
!*Staggering in X direction
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        Wt_X(1:nVar,nElemsX+ii+1,jj) = Wt_X(1:nVar,ii+1,jj) !*NB:+1 at both sides
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY    
      DO ii=1,nGhosts+1
        Wt_X(1:nVar,nElemsX+ii+1,jj) = Wt_X(1:nVar,nElemsX-ii+1,jj) !*NB:+1 at LHS
      END DO
    END DO
  CASE(3) ! Inflow
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        Wt_X(1:nVar,nElemsX+ii+1,jj) = 0.0 !*NB:+1 at LHS
      END DO
    END DO
  CASE(4) ! Outflow
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        Wt_X(1:nVar,nElemsX+ii+1,jj) = 0.0 !*NB:+1 at LHS
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        Wt_X(1:nVar,nElemsX+ii+1,jj) = Wt_X(1:nVar,nElemsX-ii+1,jj) !*NB:+1 at LHS
        Wt_X(idx_vx,nElemsX+ii+1,jj) =-Wt_X(idx_vx,nElemsX-ii+1,jj) !*NB:+1 at LHS
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition R not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!*Staggering in Y direction
!*NB: +1 in loop in Y direction
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=1,nGhosts+1
        Wt_Y(1:nVar,nElemsX+ii,jj) = Wt_Y(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction    
      DO ii=1,nGhosts+1
        Wt_Y(1:nVar,nElemsX+ii,jj) = Wt_Y(1:nVar,nElemsX-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=1,nGhosts+1
        Wt_Y(1:nVar,nElemsX+ii,jj) = 0.0
      END DO
    END DO
  CASE(4) ! Outflow
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=1,nGhosts+1
        Wt_Y(1:nVar,nElemsX+ii,jj) = 0.0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=1,nGhosts+1
        Wt_Y(1:nVar,nElemsX+ii,jj) = Wt_Y(1:nVar,nElemsX-ii+1,jj)
        Wt_Y(idx_vx,nElemsX+ii,jj) =-Wt_Y(idx_vx,nElemsX-ii+1,jj)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition R not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!------------------------------!
! Top Boundary Conditions      !
!------------------------------!
!*Staggering in X direction
!*NB: +1 in loop in X direction
SELECT CASE(BoundaryConditionsType(3))
  CASE(1) ! Periodic
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=1,nGhosts+1
        Wt_X(1:nVar,ii,nElemsY+jj) = Wt_X(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=1,nGhosts+1
        Wt_X(1:nVar,ii,nElemsY+jj) = Wt_X(1:nVar,ii,nElemsY-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=1,nGhosts+1
        Wt_X(1:nVar,ii,nElemsY+jj) = 0.0
      END DO
    END DO
  CASE(4) ! Outflow
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=1,nGhosts+1
        Wt_X(1:nVar,ii,nElemsY+jj) = 0.0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=1,nGhosts+1
        Wt_X(1:nVar,ii,nElemsY+jj) = Wt_X(1:nVar,ii,nElemsY-jj+1)
        Wt_X(idx_vy,ii,nElemsY+jj) =-Wt_X(idx_vy,ii,nElemsY-jj+1)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition T not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!*Staggering in Y direction
SELECT CASE(BoundaryConditionsType(3))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        Wt_Y(1:nVar,ii,nElemsY+jj+1) = Wt_Y(1:nVar,ii,jj+1) !*NB:+1 at both sides
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        Wt_Y(1:nVar,ii,nElemsY+jj+1) = Wt_Y(1:nVar,ii,nElemsY-jj+1) !*NB:+1 at LHS
      END DO
    END DO
  CASE(3) ! Inflow
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        Wt_Y(1:nVar,ii,nElemsY+jj+1) = 0.0 !*NB:+1 at LHS
      END DO
    END DO
  CASE(4) ! Outflow
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        Wt_Y(1:nVar,ii,nElemsY+jj+1) = 0.0 !*NB:+1 at LHS
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        Wt_Y(1:nVar,ii,nElemsY+jj+1) = Wt_Y(1:nVar,ii,nElemsY-jj+1) !"NB:+1 at LHS
        Wt_Y(idx_vy,ii,nElemsY+jj+1) =-Wt_Y(idx_vy,ii,nElemsY-jj+1) !"NB:+1 at LHS
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition T not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!------------------------------!
! Bottom Boundary Conditions   !
!------------------------------!
!*Staggering in X direction
!*NB: +1 in loop in X direction
SELECT CASE(BoundaryConditionsType(1))
  CASE(1) ! Periodic
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=0,nGhosts
        Wt_X(1:nVar,ii,-nGhosts+jj) = Wt_X(1:nVar,ii,nElemsY-nGhosts+jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=0,nGhosts
        Wt_X(1:nVar,ii,-nGhosts+jj) = Wt_X(1:nVar,ii,nGhosts-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=0,nGhosts
        Wt_X(1:nVar,ii,-nGhosts+jj) = 0.0
      END DO
    END DO
  CASE(4) ! Outflow
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=0,nGhosts
        Wt_X(1:nVar,ii,-nGhosts+jj) = 0.0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=0,nGhosts
        Wt_X(1:nVar,ii,-nGhosts+jj) = Wt_X(1:nVar,ii,nGhosts-jj+1)
        Wt_X(idx_vy,ii,-nGhosts+jj) =-Wt_X(idx_vy,ii,nGhosts-jj+1)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition B not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!*Staggering in Y direction
SELECT CASE(BoundaryConditionsType(1))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        Wt_Y(1:nVar,ii,-nGhosts+jj) = Wt_Y(1:nVar,ii,nElemsY-nGhosts+jj) !*NB: NO CHANGE IN INDICES
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        Wt_Y(1:nVar,ii,-nGhosts+jj) = Wt_Y(1:nVar,ii,nGhosts-jj+1+1) !*NB:+1 at the RHS
      END DO
    END DO
  CASE(3) ! Inflow
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        Wt_Y(1:nVar,ii,-nGhosts+jj) = 0.0
      END DO
    END DO
  CASE(4) ! Outflow
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        Wt_Y(1:nVar,ii,-nGhosts+jj) = 0.0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        Wt_Y(1:nVar,ii,-nGhosts+jj) = Wt_Y(1:nVar,ii,nGhosts-jj+1+1) !*NB: +1 at the RHS
        Wt_Y(idx_vy,ii,-nGhosts+jj) =-Wt_Y(idx_vy,ii,nGhosts-jj+1+1) !*NB: +1 at the RHS
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition B not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!------------------------------!
! Left Corners Boundary Conditions!
!------------------------------!
!*Staggering in X direction
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    !*Bottom left
    DO jj=-nGhosts,0
      DO ii=0,nGhosts
        Wt_X(1:nVar,-nGhosts+ii,jj) = Wt_X(1:nVar,nElemsX-nGhosts+ii,jj) !*NB: NO CHANGE IN INDICES
      END DO
    END DO

    !*Top left
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=0,nGhosts
        Wt_X(1:nVar,-nGhosts+ii,jj) = Wt_X(1:nVar,nElemsX-nGhosts+ii,jj) !*NB: NO CHANGE IN INDICES
      END DO
    END DO
END SELECT

!*Staggering in Y direction
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    !*Bottom left
    DO jj=-nGhosts,0
      DO ii=0,nGhosts
        Wt_Y(1:nVar,-nGhosts+ii,jj) = Wt_Y(1:nVar,nElemsX-nGhosts+ii,jj) !*NB: NO CHANGE IN INDICES
      END DO
    END DO

    !*Top left
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=0,nGhosts
        Wt_Y(1:nVar,-nGhosts+ii,jj+1) = Wt_Y(1:nVar,nElemsX-nGhosts+ii,jj+1) !*NB:+1 at both sides
      END DO
    END DO
END SELECT



!------------------------------!
! Right Corners Boundary Conditions    !
!------------------------------!
!*Staggering in X direction
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    !*Bottom right
    DO jj=-nGhosts,0
      DO ii=1,nGhosts+1
        Wt_X(1:nVar,nElemsX+ii+1,jj) = Wt_X(1:nVar,ii+1,jj) !*NB:+1 at both sides
      END DO
    END DO

    !*Top right
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=1,nGhosts+1
        Wt_X(1:nVar,nElemsX+ii+1,jj) = Wt_X(1:nVar,ii+1,jj) !*NB:+1 at both sides
      END DO
    END DO
END SELECT

!*Staggering in Y direction
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    !*Bottom right
    DO jj=-nGhosts,0
      DO ii=1,nGhosts+1
        Wt_Y(1:nVar,nElemsX+ii,jj) = Wt_Y(1:nVar,ii,jj) !*NB: NO CHANGE IN INDICES
      END DO
    END DO

    !*Top right
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=1,nGhosts+1
        Wt_Y(1:nVar,nElemsX+ii,jj+1) = Wt_Y(1:nVar,ii,jj+1) !*NB:+1 at both sides
      END DO
    END DO
END SELECT



!-------------------------------------------------------------------------------!
END SUBROUTINE Impose_BC_on_Wt
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef CENTEREDPRIMITIVE
SUBROUTINE Impose_BC_on_WCt()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsType
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
INTEGER            :: idx_vx, idx_vy
REAL               :: x0, xc, xt
REAL               :: Prim_in(1:nVar), Prim_out(1:nVar)
REAL               :: Cons_in(1:nVar), Cons_out(1:nVar)
REAL               :: ConsRefState1(1:nVar), ConsRefState2(1:nVar)
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

idx_vx = 2
idx_vy = 3

!------------------------------!
! Left Boundary Conditions     !
!------------------------------!
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        WCt(1:nVar,-nGhosts+ii,jj) = WCt(1:nVar,nElemsX-nGhosts+ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        ! PRINT*, "CORRECT BCs"
        WCt(1:nVar,-nGhosts+ii,jj) = WCt(1:nVar,nGhosts-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        WCt(1:nVar,-nGhosts+ii,jj) = 0.0
      END DO
    END DO
  CASE(4) ! Outflow
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        WCt(1:nVar,-nGhosts+ii,jj) = 0.0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        WCt(1:nVar,-nGhosts+ii,jj) = WCt(1:nVar,nGhosts-ii+1,jj)
        WCt(idx_vx,-nGhosts+ii,jj) =-WCt(idx_vx,nGhosts-ii+1,jj)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition L not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT



!------------------------------!
! Right Boundary Conditions    !
!------------------------------!
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        WCt(1:nVar,nElemsX+ii,jj) = WCt(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY    
      DO ii=1,nGhosts+1
        WCt(1:nVar,nElemsX+ii,jj) = WCt(1:nVar,nElemsX-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        WCt(1:nVar,nElemsX+ii,jj) = 0.0
      END DO
    END DO
  CASE(4) ! Outflow
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        WCt(1:nVar,nElemsX+ii,jj) = 0.0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        WCt(1:nVar,nElemsX+ii,jj) = WCt(1:nVar,nElemsX-ii+1,jj)
        WCt(idx_vx,nElemsX+ii,jj) =-WCt(idx_vx,nElemsX-ii+1,jj)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition R not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!------------------------------!
! Top Boundary Conditions      !
!------------------------------!
SELECT CASE(BoundaryConditionsType(3))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        WCt(1:nVar,ii,nElemsY+jj) = WCt(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        WCt(1:nVar,ii,nElemsY+jj) = WCt(1:nVar,ii,nElemsY-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        WCt(1:nVar,ii,nElemsY+jj) = 0.0
      END DO
    END DO
  CASE(4) ! Outflow
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        WCt(1:nVar,ii,nElemsY+jj) = 0.0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        WCt(1:nVar,ii,nElemsY+jj) = WCt(1:nVar,ii,nElemsY-jj+1)
        WCt(idx_vy,ii,nElemsY+jj) =-WCt(idx_vy,ii,nElemsY-jj+1)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition T not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT




!------------------------------!
! Bottom Boundary Conditions   !
!------------------------------!
SELECT CASE(BoundaryConditionsType(1))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        WCt(1:nVar,ii,-nGhosts+jj) = WCt(1:nVar,ii,nElemsY-nGhosts+jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        WCt(1:nVar,ii,-nGhosts+jj) = WCt(1:nVar,ii,nGhosts-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        WCt(1:nVar,ii,-nGhosts+jj) = 0.0
      END DO
    END DO
  CASE(4) ! Outflow
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        WCt(1:nVar,ii,-nGhosts+jj) = 0.0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        WCt(1:nVar,ii,-nGhosts+jj) = WCt(1:nVar,ii,nGhosts-jj+1)
        WCt(idx_vy,ii,-nGhosts+jj) =-WCt(idx_vy,ii,nGhosts-jj+1)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition B not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT



!------------------------------!
! Left Corners Boundary Conditions!
!------------------------------!
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    !*Bottom left
    DO jj=-nGhosts,0
      DO ii=0,nGhosts
        WCt(1:nVar,-nGhosts+ii,jj) = WCt(1:nVar,nElemsX-nGhosts+ii,jj)
      END DO
    END DO

    !*Top left
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=0,nGhosts
        WCt(1:nVar,-nGhosts+ii,jj) = WCt(1:nVar,nElemsX-nGhosts+ii,jj)
      END DO
    END DO
END SELECT



!------------------------------!
! Right Corners Boundary Conditions    !
!------------------------------!
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    !*Bottom right
    DO jj=-nGhosts,0
      DO ii=1,nGhosts+1
        WCt(1:nVar,nElemsX+ii,jj) = WCt(1:nVar,ii,jj)
      END DO
    END DO

    !*Top right
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=1,nGhosts+1
        WCt(1:nVar,nElemsX+ii,jj) = WCt(1:nVar,ii,jj)
      END DO
    END DO
END SELECT



!-------------------------------------------------------------------------------!
END SUBROUTINE Impose_BC_on_WCt
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE Impose_BC_on_W(WX,WY,t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsType
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: V
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(INOUT)    :: WX(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
REAL,INTENT(INOUT)    :: WY(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
REAL,INTENT(IN)    :: t

!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
INTEGER            :: idx_vx, idx_vy
REAL               :: x0, xc, xt
REAL               :: Prim_in(1:nVar), Prim_out(1:nVar)
REAL               :: Cons_in(1:nVar), Cons_out(1:nVar)
REAL               :: ConsRefState1(1:nVar), ConsRefState2(1:nVar)
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

idx_vx = 2
idx_vy = 3

!------------------------------!
! Left Boundary Conditions     !
!------------------------------!
!*Staggering in X direction
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=0,nGhosts 
        WX(1:nVar,-nGhosts+ii,jj) = WX(1:nVar,nElemsX-nGhosts+ii,jj) !*NB: NO CHANGE IN INDICES
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        WX(1:nVar,-nGhosts+ii,jj) = WX(1:nVar,nGhosts-ii+1+1,jj) !*NB: +1 at the RHS
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState4(1:nVar)
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        WX(1:nVar,-nGhosts+ii,jj) = Prim_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState4(1:nVar)
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        WX(1:nVar,-nGhosts+ii,jj) = Prim_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        WX(1:nVar,-nGhosts+ii,jj) = WX(1:nVar,nGhosts-ii+1+1,jj) !*NB:+1 at the RHS
        WX(idx_vx,-nGhosts+ii,jj) =-WX(idx_vx,nGhosts-ii+1+1,jj) !*NB:+1 at the RHS
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition L not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!*Staggering in Y direction
!*NB: +1 in loop in Y direction
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=0,nGhosts
        WY(1:nVar,-nGhosts+ii,jj) = WY(1:nVar,nElemsX-nGhosts+ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=0,nGhosts
        WY(1:nVar,-nGhosts+ii,jj) = WY(1:nVar,nGhosts-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState4(1:nVar)
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=0,nGhosts
        WY(1:nVar,-nGhosts+ii,jj) = Prim_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState4(1:nVar)
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=0,nGhosts
        WY(1:nVar,-nGhosts+ii,jj) = Prim_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=0,nGhosts
        WY(1:nVar,-nGhosts+ii,jj) = WY(1:nVar,nGhosts-ii+1,jj)
        WY(idx_vx,-nGhosts+ii,jj) =-WY(idx_vx,nGhosts-ii+1,jj)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition L not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT



!------------------------------!
! Right Boundary Conditions    !
!------------------------------!
!*Staggering in X direction
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        WX(1:nVar,nElemsX+ii+1,jj) = WX(1:nVar,ii+1,jj) !*NB:+1 at both sides
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY    
      DO ii=1,nGhosts+1
        WX(1:nVar,nElemsX+ii+1,jj) = WX(1:nVar,nElemsX-ii+1,jj) !*NB:+1 at LHS
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState2(1:nVar)
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        WX(1:nVar,nElemsX+ii+1,jj) = Prim_in(1:nVar) !*NB:+1 at LHS
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState2(1:nVar)
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        WX(1:nVar,nElemsX+ii+1,jj) = Prim_out(1:nVar) !*NB:+1 at LHS
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        WX(1:nVar,nElemsX+ii+1,jj) = WX(1:nVar,nElemsX-ii+1,jj) !*NB:+1 at LHS
        WX(idx_vx,nElemsX+ii+1,jj) =-WX(idx_vx,nElemsX-ii+1,jj) !*NB:+1 at LHS
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition R not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!*Staggering in Y direction
!*NB: +1 in loop in Y direction
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=1,nGhosts+1
        WY(1:nVar,nElemsX+ii,jj) = WY(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction    
      DO ii=1,nGhosts+1
        WY(1:nVar,nElemsX+ii,jj) = WY(1:nVar,nElemsX-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState2(1:nVar)
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=1,nGhosts+1
        WY(1:nVar,nElemsX+ii,jj) = Prim_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState2(1:nVar)
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=1,nGhosts+1
        WY(1:nVar,nElemsX+ii,jj) = Prim_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=1,nGhosts+1
        WY(1:nVar,nElemsX+ii,jj) = WY(1:nVar,nElemsX-ii+1,jj)
        WY(idx_vx,nElemsX+ii,jj) =-WY(idx_vx,nElemsX-ii+1,jj)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition R not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!------------------------------!
! Top Boundary Conditions      !
!------------------------------!
!*Staggering in X direction
!*NB: +1 in loop in X direction
SELECT CASE(BoundaryConditionsType(3))
  CASE(1) ! Periodic
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=1,nGhosts+1
        WX(1:nVar,ii,nElemsY+jj) = WX(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=1,nGhosts+1
        WX(1:nVar,ii,nElemsY+jj) = WX(1:nVar,ii,nElemsY-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState3(1:nVar)
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=1,nGhosts+1
        WX(1:nVar,ii,nElemsY+jj) = Prim_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState3(1:nVar)
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=1,nGhosts+1
        WX(1:nVar,ii,nElemsY+jj) = Prim_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=1,nGhosts+1
        WX(1:nVar,ii,nElemsY+jj) = WX(1:nVar,ii,nElemsY-jj+1)
        WX(idx_vy,ii,nElemsY+jj) =-WX(idx_vy,ii,nElemsY-jj+1)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition T not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!*Staggering in Y direction
SELECT CASE(BoundaryConditionsType(3))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        WY(1:nVar,ii,nElemsY+jj+1) = WY(1:nVar,ii,jj+1) !*NB:+1 at both sides
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        WY(1:nVar,ii,nElemsY+jj+1) = WY(1:nVar,ii,nElemsY-jj+1) !*NB:+1 at LHS
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState3(1:nVar)
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        WY(1:nVar,ii,nElemsY+jj+1) = Prim_in(1:nVar) !*NB:+1 at LHS
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState3(1:nVar)
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        WY(1:nVar,ii,nElemsY+jj+1) = Prim_out(1:nVar) !*NB:+1 at LHS
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        WY(1:nVar,ii,nElemsY+jj+1) = WY(1:nVar,ii,nElemsY-jj+1) !"NB:+1 at LHS
        WY(idx_vy,ii,nElemsY+jj+1) =-WY(idx_vy,ii,nElemsY-jj+1) !"NB:+1 at LHS
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition T not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!------------------------------!
! Bottom Boundary Conditions   !
!------------------------------!
!*Staggering in X direction
!*NB: +1 in loop in X direction
SELECT CASE(BoundaryConditionsType(1))
  CASE(1) ! Periodic
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=0,nGhosts
        WX(1:nVar,ii,-nGhosts+jj) = WX(1:nVar,ii,nElemsY-nGhosts+jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=0,nGhosts
        WX(1:nVar,ii,-nGhosts+jj) = WX(1:nVar,ii,nGhosts-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState1(1:nVar)
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=0,nGhosts
        WX(1:nVar,ii,-nGhosts+jj) = Prim_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState1(1:nVar)
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=0,nGhosts
        WX(1:nVar,ii,-nGhosts+jj) = Prim_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=0,nGhosts
        WX(1:nVar,ii,-nGhosts+jj) = WX(1:nVar,ii,nGhosts-jj+1)
        WX(idx_vy,ii,-nGhosts+jj) =-WX(idx_vy,ii,nGhosts-jj+1)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition B not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!*Staggering in Y direction
SELECT CASE(BoundaryConditionsType(1))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        WY(1:nVar,ii,-nGhosts+jj) = WY(1:nVar,ii,nElemsY-nGhosts+jj) !*NB: NO CHANGE IN INDICES
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        WY(1:nVar,ii,-nGhosts+jj) = WY(1:nVar,ii,nGhosts-jj+1+1) !*NB:+1 at the RHS
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState1(1:nVar)
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        WY(1:nVar,ii,-nGhosts+jj) = Prim_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState1(1:nVar)
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        WY(1:nVar,ii,-nGhosts+jj) = Prim_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        WY(1:nVar,ii,-nGhosts+jj) = WY(1:nVar,ii,nGhosts-jj+1+1) !*NB: +1 at the RHS
        WY(idx_vy,ii,-nGhosts+jj) =-WY(idx_vy,ii,nGhosts-jj+1+1) !*NB: +1 at the RHS
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition B not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!------------------------------!
! Left Corners Boundary Conditions!
!------------------------------!
!*Staggering in X direction
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    !*Bottom left
    DO jj=-nGhosts,0
      DO ii=0,nGhosts
        WX(1:nVar,-nGhosts+ii,jj) = WX(1:nVar,nElemsX-nGhosts+ii,jj) !*NB: NO CHANGE IN INDICES
      END DO
    END DO

    !*Top left
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=0,nGhosts
        WX(1:nVar,-nGhosts+ii,jj) = WX(1:nVar,nElemsX-nGhosts+ii,jj) !*NB: NO CHANGE IN INDICES
      END DO
    END DO
END SELECT

!*Staggering in Y direction
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    !*Bottom left
    DO jj=-nGhosts,0
      DO ii=0,nGhosts
        WY(1:nVar,-nGhosts+ii,jj) = WY(1:nVar,nElemsX-nGhosts+ii,jj) !*NB: NO CHANGE IN INDICES
      END DO
    END DO

    !*Top left
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=0,nGhosts
        WY(1:nVar,-nGhosts+ii,jj+1) = WY(1:nVar,nElemsX-nGhosts+ii,jj+1) !*NB:+1 at both sides
      END DO
    END DO
END SELECT




!------------------------------!
! Right Corners Boundary Conditions    !
!------------------------------!
!*Staggering in X direction
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    !*Bottom right
    DO jj=-nGhosts,0
      DO ii=1,nGhosts+1
        WX(1:nVar,nElemsX+ii+1,jj) = WX(1:nVar,ii+1,jj) !*NB:+1 at both sides
      END DO
    END DO

    !*Top right
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=1,nGhosts+1
        WX(1:nVar,nElemsX+ii+1,jj) = WX(1:nVar,ii+1,jj) !*NB:+1 at both sides
      END DO
    END DO
END SELECT

!*Staggering in Y direction
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    !*Bottom right
    DO jj=-nGhosts,0
      DO ii=1,nGhosts+1
        WY(1:nVar,nElemsX+ii,jj) = WY(1:nVar,ii,jj) !*NB: NO CHANGE IN INDICES
      END DO
    END DO

    !*Top right
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=1,nGhosts+1
        WY(1:nVar,nElemsX+ii,jj+1) = WY(1:nVar,ii,jj+1) !*NB:+1 at both sides
      END DO
    END DO
END SELECT



!-------------------------------------------------------------------------------!
END SUBROUTINE Impose_BC_on_W
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef CENTEREDPRIMITIVE
SUBROUTINE Impose_BC_on_WC(WC,t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsType
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(INOUT)    :: WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
REAL,INTENT(IN)       :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
INTEGER            :: idx_vx, idx_vy
REAL               :: x0, xc, xt
REAL               :: Prim_in(1:nVar), Prim_out(1:nVar)
REAL               :: Cons_in(1:nVar), Cons_out(1:nVar)
REAL               :: ConsRefState1(1:nVar), ConsRefState2(1:nVar)
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

idx_vx = 2
idx_vy = 3

!------------------------------!
! Left Boundary Conditions     !
!------------------------------!
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        WC(1:nVar,-nGhosts+ii,jj) = WC(1:nVar,nElemsX-nGhosts+ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        WC(1:nVar,-nGhosts+ii,jj) = WC(1:nVar,nGhosts-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState4(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        WC(1:nVar,-nGhosts+ii,jj) = Prim_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState4(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        WC(1:nVar,-nGhosts+ii,jj) = Prim_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        WC(1:nVar,-nGhosts+ii,jj) = WC(1:nVar,nGhosts-ii+1,jj)
        WC(idx_vx,-nGhosts+ii,jj) =-WC(idx_vx,nGhosts-ii+1,jj)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition L not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT



!------------------------------!
! Right Boundary Conditions    !
!------------------------------!
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        WC(1:nVar,nElemsX+ii,jj) = WC(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY    
      DO ii=1,nGhosts+1
        WC(1:nVar,nElemsX+ii,jj) = WC(1:nVar,nElemsX-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState2(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        WC(1:nVar,nElemsX+ii,jj) = Prim_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState2(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        WC(1:nVar,nElemsX+ii,jj) = Prim_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        WC(1:nVar,nElemsX+ii,jj) = WC(1:nVar,nElemsX-ii+1,jj)
        WC(idx_vx,nElemsX+ii,jj) =-WC(idx_vx,nElemsX-ii+1,jj)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition R not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!------------------------------!
! Top Boundary Conditions      !
!------------------------------!
SELECT CASE(BoundaryConditionsType(3))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        WC(1:nVar,ii,nElemsY+jj) = WC(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        WC(1:nVar,ii,nElemsY+jj) = WC(1:nVar,ii,nElemsY-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState3(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        WC(1:nVar,ii,nElemsY+jj) = Prim_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState3(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        WC(1:nVar,ii,nElemsY+jj) = Prim_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        WC(1:nVar,ii,nElemsY+jj) = WC(1:nVar,ii,nElemsY-jj+1)
        WC(idx_vy,ii,nElemsY+jj) =-WC(idx_vy,ii,nElemsY-jj+1)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition T not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT



!------------------------------!
! Bottom Boundary Conditions   !
!------------------------------!
SELECT CASE(BoundaryConditionsType(1))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        WC(1:nVar,ii,-nGhosts+jj) = WC(1:nVar,ii,nElemsY-nGhosts+jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        WC(1:nVar,ii,-nGhosts+jj) = WC(1:nVar,ii,nGhosts-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState1(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        WC(1:nVar,ii,-nGhosts+jj) = Prim_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState1(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        WC(1:nVar,ii,-nGhosts+jj) = Prim_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        WC(1:nVar,ii,-nGhosts+jj) = WC(1:nVar,ii,nGhosts-jj+1)
        WC(idx_vy,ii,-nGhosts+jj) =-WC(idx_vy,ii,nGhosts-jj+1)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition B not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT



!------------------------------!
! Left Corners Boundary Conditions!
!------------------------------!
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    !*Bottom left
    DO jj=-nGhosts,0
      DO ii=0,nGhosts
        WC(1:nVar,-nGhosts+ii,jj) = WC(1:nVar,nElemsX-nGhosts+ii,jj)
      END DO
    END DO

    !*Top left
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=0,nGhosts
        WC(1:nVar,-nGhosts+ii,jj) = WC(1:nVar,nElemsX-nGhosts+ii,jj)
      END DO
    END DO
END SELECT



!------------------------------!
! Right Corners Boundary Conditions    !
!------------------------------!
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    !*Bottom right
    DO jj=-nGhosts,0
      DO ii=1,nGhosts+1
        WC(1:nVar,nElemsX+ii,jj) = WC(1:nVar,ii,jj)
      END DO
    END DO

    !*Top right
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=1,nGhosts+1
        WC(1:nVar,nElemsX+ii,jj) = WC(1:nVar,ii,jj)
      END DO
    END DO
END SELECT



!-------------------------------------------------------------------------------!
END SUBROUTINE Impose_BC_on_WC
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX)
#ifdef MULTIFLUID
REAL FUNCTION Get_Gamma(level_set_variable)
  USE MOD_FiniteVolume2D_vars,ONLY: Gmm1
  USE MOD_FiniteVolume2D_vars,ONLY: Gmm2
  IMPLICIT NONE
  REAL, INTENT(IN) :: level_set_variable
  REAL, PARAMETER  :: smoothing_level_set_variable=1e-3

  IF (level_set_variable .GE. 0.0) THEN
    Get_Gamma=Gmm1
  ELSE
    Get_Gamma=Gmm2
  END IF

  ! IF (level_set_variable .GE. smoothing_level_set_variable) THEN
  !   Get_Gamma=Gmm1
  ! ELSE IF (level_set_variable .LE. -smoothing_level_set_variable) THEN
  !   Get_Gamma=Gmm2
  ! ELSE
  !   Get_Gamma=Gmm2+(Gmm1-Gmm2)/(2.0*smoothing_level_set_variable)*(level_set_variable+smoothing_level_set_variable)
  ! END IF

END FUNCTION Get_Gamma
#endif
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX)
#ifdef MULTIFLUID
REAL FUNCTION Get_Pressure_Infinity(level_set_variable)
  USE MOD_FiniteVolume2D_vars,ONLY: pinf1
  USE MOD_FiniteVolume2D_vars,ONLY: pinf2
  IMPLICIT NONE
  REAL, INTENT(IN) :: level_set_variable
  REAL, PARAMETER  :: smoothing_level_set_variable=1e-3

  IF (level_set_variable .GE. 0.0) THEN
    Get_Pressure_Infinity=pinf1
  ELSE
    Get_Pressure_Infinity=pinf2
  END IF

  ! IF (level_set_variable .GE. smoothing_level_set_variable) THEN
  !   Get_Pressure_Infinity=pinf1
  ! ELSE IF (level_set_variable .LE. -smoothing_level_set_variable) THEN
  !   Get_Pressure_Infinity=pinf2
  ! ELSE
  !   Get_Pressure_Infinity=pinf2+(pinf1-pinf2)/(2.0*smoothing_level_set_variable)*(level_set_variable+smoothing_level_set_variable)
  ! END IF

END FUNCTION Get_Pressure_Infinity
#endif
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX)
#ifdef MULTIFLUID
SUBROUTINE Impose_BC_on_Troubled_Cell_U()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsType
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_U
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
REAL               :: x0, xc, xt
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!


!------------------------------!
! Left Boundary Conditions     !
!------------------------------!
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        Troubled_Cell_U(-nGhosts+ii,jj) = Troubled_Cell_U(nElemsX-nGhosts+ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        Troubled_Cell_U(-nGhosts+ii,jj) = 0!*Troubled_Cell_U(nGhosts-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        Troubled_Cell_U(-nGhosts+ii,jj) = 0
      END DO
    END DO
  CASE(4) ! Outflow
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        Troubled_Cell_U(-nGhosts+ii,jj) = 0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        Troubled_Cell_U(-nGhosts+ii,jj) = 0!*Troubled_Cell_U(nGhosts-ii+1,jj)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition L not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!------------------------------!
! Right Boundary Conditions    !
!------------------------------!
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        Troubled_Cell_U(nElemsX+ii,jj) = Troubled_Cell_U(ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY    
      DO ii=1,nGhosts+1
        Troubled_Cell_U(nElemsX+ii,jj) = 0!*Troubled_Cell_U(nElemsX-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        Troubled_Cell_U(nElemsX+ii,jj) = 0
      END DO
    END DO
  CASE(4) ! Outflow
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        Troubled_Cell_U(nElemsX+ii,jj) = 0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        Troubled_Cell_U(nElemsX+ii,jj) = 0!*Troubled_Cell_U(nElemsX-ii+1,jj)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition R not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT



!------------------------------!
! Top Boundary Conditions      !
!------------------------------!
SELECT CASE(BoundaryConditionsType(3))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        Troubled_Cell_U(ii,nElemsY+jj) = Troubled_Cell_U(ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        Troubled_Cell_U(ii,nElemsY+jj) = 0!*Troubled_Cell_U(ii,nElemsY-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        Troubled_Cell_U(ii,nElemsY+jj) = 0
      END DO
    END DO
  CASE(4) ! Outflow
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        Troubled_Cell_U(ii,nElemsY+jj) = 0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        Troubled_Cell_U(ii,nElemsY+jj) = 0!*Troubled_Cell_U(ii,nElemsY-jj+1)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition T not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT



!------------------------------!
! Bottom Boundary Conditions   !
!------------------------------!
SELECT CASE(BoundaryConditionsType(1))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        Troubled_Cell_U(ii,-nGhosts+jj) = Troubled_Cell_U(ii,nElemsY-nGhosts+jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        Troubled_Cell_U(ii,-nGhosts+jj) = 0!*Troubled_Cell_U(ii,nGhosts-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        Troubled_Cell_U(ii,-nGhosts+jj) = 0
      END DO
    END DO
  CASE(4) ! Outflow
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        Troubled_Cell_U(ii,-nGhosts+jj) = 0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        Troubled_Cell_U(ii,-nGhosts+jj) = 0!*Troubled_Cell_U(ii,nGhosts-jj+1)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition B not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT



!------------------------------!
! Left Corners Boundary Conditions!
!------------------------------!
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    !*Bottom left
    DO jj=-nGhosts,0
      DO ii=0,nGhosts
        Troubled_Cell_U(-nGhosts+ii,jj) = Troubled_Cell_U(nElemsX-nGhosts+ii,jj)
      END DO
    END DO

    !*Top left
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=0,nGhosts
        Troubled_Cell_U(-nGhosts+ii,jj) = Troubled_Cell_U(nElemsX-nGhosts+ii,jj)
      END DO
    END DO
END SELECT



!------------------------------!
! Right Corners Boundary Conditions    !
!------------------------------!
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    !*Bottom right
    DO jj=-nGhosts,0
      DO ii=1,nGhosts+1
        Troubled_Cell_U(nElemsX+ii,jj) = Troubled_Cell_U(ii,jj)
      END DO
    END DO

    !*Top right
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=1,nGhosts+1
        Troubled_Cell_U(nElemsX+ii,jj) = Troubled_Cell_U(ii,jj)
      END DO
    END DO
END SELECT



!-------------------------------------------------------------------------------!
END SUBROUTINE Impose_BC_on_Troubled_Cell_U
#endif
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(CENTEREDPRIMITIVE) || defined(PATHCONSERVATIVESHOCKDETECTION)
SUBROUTINE Impose_BC_on_Troubled_Cell_INPUT(Troubled_Cell)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsType
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
INTEGER, INTENT(INOUT) :: Troubled_Cell(-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) 
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
REAL               :: x0, xc, xt
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!


!------------------------------!
! Left Boundary Conditions     !
!------------------------------!
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        Troubled_Cell(-nGhosts+ii,jj) = Troubled_Cell(nElemsX-nGhosts+ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        Troubled_Cell(-nGhosts+ii,jj) = 0
      END DO
    END DO
  CASE(3) ! Inflow
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        Troubled_Cell(-nGhosts+ii,jj) = 0
      END DO
    END DO
  CASE(4) ! Outflow
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        Troubled_Cell(-nGhosts+ii,jj) = 0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        Troubled_Cell(-nGhosts+ii,jj) = 0
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition L not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!------------------------------!
! Right Boundary Conditions    !
!------------------------------!
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        Troubled_Cell(nElemsX+ii,jj) = Troubled_Cell(ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY    
      DO ii=1,nGhosts+1
        Troubled_Cell(nElemsX+ii,jj) = 0
      END DO
    END DO
  CASE(3) ! Inflow
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        Troubled_Cell(nElemsX+ii,jj) = 0
      END DO
    END DO
  CASE(4) ! Outflow
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        Troubled_Cell(nElemsX+ii,jj) = 0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        Troubled_Cell(nElemsX+ii,jj) = 0
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition R not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT



!------------------------------!
! Top Boundary Conditions      !
!------------------------------!
SELECT CASE(BoundaryConditionsType(3))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        Troubled_Cell(ii,nElemsY+jj) = Troubled_Cell(ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        Troubled_Cell(ii,nElemsY+jj) = 0
      END DO
    END DO
  CASE(3) ! Inflow
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        Troubled_Cell(ii,nElemsY+jj) = 0
      END DO
    END DO
  CASE(4) ! Outflow
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        Troubled_Cell(ii,nElemsY+jj) = 0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        Troubled_Cell(ii,nElemsY+jj) = 0
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition T not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT



!------------------------------!
! Bottom Boundary Conditions   !
!------------------------------!
SELECT CASE(BoundaryConditionsType(1))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        Troubled_Cell(ii,-nGhosts+jj) = Troubled_Cell(ii,nElemsY-nGhosts+jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        Troubled_Cell(ii,-nGhosts+jj) = 0
      END DO
    END DO
  CASE(3) ! Inflow
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        Troubled_Cell(ii,-nGhosts+jj) = 0
      END DO
    END DO
  CASE(4) ! Outflow
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        Troubled_Cell(ii,-nGhosts+jj) = 0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        Troubled_Cell(ii,-nGhosts+jj) = 0
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition B not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT



!------------------------------!
! Left Corners Boundary Conditions!
!------------------------------!
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    !*Bottom left
    DO jj=-nGhosts,0
      DO ii=0,nGhosts
        Troubled_Cell(-nGhosts+ii,jj) = Troubled_Cell(nElemsX-nGhosts+ii,jj)
      END DO
    END DO

    !*Top left
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=0,nGhosts
        Troubled_Cell(-nGhosts+ii,jj) = Troubled_Cell(nElemsX-nGhosts+ii,jj)
      END DO
    END DO
END SELECT



!------------------------------!
! Right Corners Boundary Conditions    !
!------------------------------!
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    !*Bottom right
    DO jj=-nGhosts,0
      DO ii=1,nGhosts+1
        Troubled_Cell(nElemsX+ii,jj) = Troubled_Cell(ii,jj)
      END DO
    END DO

    !*Top right
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=1,nGhosts+1
        Troubled_Cell(nElemsX+ii,jj) = Troubled_Cell(ii,jj)
      END DO
    END DO
END SELECT



!-------------------------------------------------------------------------------!
END SUBROUTINE Impose_BC_on_Troubled_Cell_INPUT
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX)
#ifdef MULTIFLUID
SUBROUTINE Impose_BC_on_Troubled_Cell_W()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsType
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_W_X
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_W_Y
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
INTEGER            :: idx_vx, idx_vy
REAL               :: x0, xc, xt
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

idx_vx = 2
idx_vy = 3

!------------------------------!
! Left Boundary Conditions     !
!------------------------------!
!*Staggering in X direction
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=0,nGhosts 
        Troubled_Cell_W_X(-nGhosts+ii,jj) = Troubled_Cell_W_X(nElemsX-nGhosts+ii,jj) !*NB: NO CHANGE IN INDICES
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        Troubled_Cell_W_X(-nGhosts+ii,jj) = 0
      END DO
    END DO
  CASE(3) ! Inflow
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        Troubled_Cell_W_X(-nGhosts+ii,jj) = 0
      END DO
    END DO
  CASE(4) ! Outflow
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        Troubled_Cell_W_X(-nGhosts+ii,jj) = 0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        Troubled_Cell_W_X(-nGhosts+ii,jj) = 0
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition L not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!*Staggering in Y direction
!*NB: +1 in loop in Y direction
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=0,nGhosts
        Troubled_Cell_W_Y(-nGhosts+ii,jj) = Troubled_Cell_W_Y(nElemsX-nGhosts+ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=0,nGhosts
        Troubled_Cell_W_Y(-nGhosts+ii,jj) = 0
      END DO
    END DO
  CASE(3) ! Inflow
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=0,nGhosts
        Troubled_Cell_W_Y(-nGhosts+ii,jj) = 0
      END DO
    END DO
  CASE(4) ! Outflow
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=0,nGhosts
        Troubled_Cell_W_Y(-nGhosts+ii,jj) = 0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=0,nGhosts
        Troubled_Cell_W_Y(-nGhosts+ii,jj) = 0
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition L not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT



!------------------------------!
! Right Boundary Conditions    !
!------------------------------!
!*Staggering in X direction
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        Troubled_Cell_W_X(nElemsX+ii+1,jj) = Troubled_Cell_W_X(ii+1,jj) !*NB:+1 at both sides
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY    
      DO ii=1,nGhosts+1
        Troubled_Cell_W_X(nElemsX+ii+1,jj) = 0
      END DO
    END DO
  CASE(3) ! Inflow
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        Troubled_Cell_W_X(nElemsX+ii+1,jj) = 0
      END DO
    END DO
  CASE(4) ! Outflow
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        Troubled_Cell_W_X(nElemsX+ii+1,jj) = 0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        Troubled_Cell_W_X(nElemsX+ii+1,jj) = 0
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition R not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!*Staggering in Y direction
!*NB: +1 in loop in Y direction
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=1,nGhosts+1
        Troubled_Cell_W_Y(nElemsX+ii,jj) = Troubled_Cell_W_Y(ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction    
      DO ii=1,nGhosts+1
        Troubled_Cell_W_Y(nElemsX+ii,jj) = 0
      END DO
    END DO
  CASE(3) ! Inflow
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=1,nGhosts+1
        Troubled_Cell_W_Y(nElemsX+ii,jj) = 0
      END DO
    END DO
  CASE(4) ! Outflow
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=1,nGhosts+1
        Troubled_Cell_W_Y(nElemsX+ii,jj) = 0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY+1 !*NB: +1 in Y direction
      DO ii=1,nGhosts+1
        Troubled_Cell_W_Y(nElemsX+ii,jj) = 0
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition R not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!------------------------------!
! Top Boundary Conditions      !
!------------------------------!
!*Staggering in X direction
!*NB: +1 in loop in X direction
SELECT CASE(BoundaryConditionsType(3))
  CASE(1) ! Periodic
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=1,nGhosts+1
        Troubled_Cell_W_X(ii,nElemsY+jj) = Troubled_Cell_W_X(ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=1,nGhosts+1
        Troubled_Cell_W_X(ii,nElemsY+jj) = 0
      END DO
    END DO
  CASE(3) ! Inflow
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=1,nGhosts+1
        Troubled_Cell_W_X(ii,nElemsY+jj) = 0
      END DO
    END DO
  CASE(4) ! Outflow
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=1,nGhosts+1
        Troubled_Cell_W_X(ii,nElemsY+jj) = 0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=1,nGhosts+1
        Troubled_Cell_W_X(ii,nElemsY+jj) = 0
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition T not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!*Staggering in Y direction
SELECT CASE(BoundaryConditionsType(3))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        Troubled_Cell_W_Y(ii,nElemsY+jj+1) = Troubled_Cell_W_Y(ii,jj+1) !*NB:+1 at both sides
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        Troubled_Cell_W_Y(ii,nElemsY+jj+1) = 0
      END DO
    END DO
  CASE(3) ! Inflow
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        Troubled_Cell_W_Y(ii,nElemsY+jj+1) = 0
      END DO
    END DO
  CASE(4) ! Outflow
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        Troubled_Cell_W_Y(ii,nElemsY+jj+1) = 0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        Troubled_Cell_W_Y(ii,nElemsY+jj+1) = 0
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition T not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!------------------------------!
! Bottom Boundary Conditions   !
!------------------------------!
!*Staggering in X direction
!*NB: +1 in loop in X direction
SELECT CASE(BoundaryConditionsType(1))
  CASE(1) ! Periodic
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=0,nGhosts
        Troubled_Cell_W_X(ii,-nGhosts+jj) = Troubled_Cell_W_X(ii,nElemsY-nGhosts+jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=0,nGhosts
        Troubled_Cell_W_X(ii,-nGhosts+jj) = 0
      END DO
    END DO
  CASE(3) ! Inflow
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=0,nGhosts
        Troubled_Cell_W_X(ii,-nGhosts+jj) = 0
      END DO
    END DO
  CASE(4) ! Outflow
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=0,nGhosts
        Troubled_Cell_W_X(ii,-nGhosts+jj) = 0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX+1 !*NB: +1 in X direction
      DO jj=0,nGhosts
        Troubled_Cell_W_X(ii,-nGhosts+jj) = 0
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition B not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!*Staggering in Y direction
SELECT CASE(BoundaryConditionsType(1))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        Troubled_Cell_W_Y(ii,-nGhosts+jj) = Troubled_Cell_W_Y(ii,nElemsY-nGhosts+jj) !*NB: NO CHANGE IN INDICES
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        Troubled_Cell_W_Y(ii,-nGhosts+jj) = 0
      END DO
    END DO
  CASE(3) ! Inflow
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        Troubled_Cell_W_Y(ii,-nGhosts+jj) = 0
      END DO
    END DO
  CASE(4) ! Outflow
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        Troubled_Cell_W_Y(ii,-nGhosts+jj) = 0
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        Troubled_Cell_W_Y(ii,-nGhosts+jj) = 0
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition B not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!------------------------------!
! Left Corners Boundary Conditions!
!------------------------------!
!*Staggering in X direction
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    !*Bottom left
    DO jj=-nGhosts,0
      DO ii=0,nGhosts
        Troubled_Cell_W_X(-nGhosts+ii,jj) = Troubled_Cell_W_X(nElemsX-nGhosts+ii,jj) !*NB: NO CHANGE IN INDICES
      END DO
    END DO

    !*Top left
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=0,nGhosts
        Troubled_Cell_W_X(-nGhosts+ii,jj) = 0
      END DO
    END DO
END SELECT

!*Staggering in Y direction
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    !*Bottom left
    DO jj=-nGhosts,0
      DO ii=0,nGhosts
        Troubled_Cell_W_Y(-nGhosts+ii,jj) = 0
      END DO
    END DO

    !*Top left
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=0,nGhosts
        Troubled_Cell_W_Y(-nGhosts+ii,jj+1) = 0
      END DO
    END DO
END SELECT




!------------------------------!
! Right Corners Boundary Conditions    !
!------------------------------!
!*Staggering in X direction
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    !*Bottom right
    DO jj=-nGhosts,0
      DO ii=1,nGhosts+1
        Troubled_Cell_W_X(nElemsX+ii+1,jj) = Troubled_Cell_W_X(ii+1,jj) !*NB:+1 at both sides
      END DO
    END DO

    !*Top right
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=1,nGhosts+1
        Troubled_Cell_W_X(nElemsX+ii+1,jj) = 0
      END DO
    END DO
END SELECT

!*Staggering in Y direction
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    !*Bottom right
    DO jj=-nGhosts,0
      DO ii=1,nGhosts+1
        Troubled_Cell_W_Y(nElemsX+ii,jj) = 0
      END DO
    END DO

    !*Top right
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=1,nGhosts+1
        Troubled_Cell_W_Y(nElemsX+ii,jj+1) = 0
      END DO
    END DO
END SELECT



!-------------------------------------------------------------------------------!
END SUBROUTINE Impose_BC_on_Troubled_Cell_W
#endif
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX)
#ifdef MULTIFLUID
SUBROUTINE Impose_BC_on_Fluid_Cell_U()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsType
USE MOD_FiniteVolume2D_vars,ONLY: Which_Fluid_In_Cell_U
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
INTEGER            :: fluid_status
REAL               :: x0, xc, xt
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!


!------------------------------!
! Left Boundary Conditions     !
!------------------------------!
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        Which_Fluid_In_Cell_U(-nGhosts+ii,jj) = Which_Fluid_In_Cell_U(nElemsX-nGhosts+ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        Which_Fluid_In_Cell_U(-nGhosts+ii,jj) = Which_Fluid_In_Cell_U(nGhosts-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    fluid_status=INT( SIGN(1.0,PrimRefState4(nVar)) )
    IF (fluid_status .EQ. 0) THEN
      fluid_status=1
    END IF
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        Which_Fluid_In_Cell_U(-nGhosts+ii,jj) = fluid_status
      END DO
    END DO
  CASE(4) ! Outflow
    fluid_status=INT( SIGN(1.0,PrimRefState4(nVar)) )
    IF (fluid_status .EQ. 0) THEN
      fluid_status=1
    END IF
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        Which_Fluid_In_Cell_U(-nGhosts+ii,jj) = fluid_status
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=0,nGhosts
        Which_Fluid_In_Cell_U(-nGhosts+ii,jj) = Which_Fluid_In_Cell_U(nGhosts-ii+1,jj)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition L not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!------------------------------!
! Right Boundary Conditions    !
!------------------------------!
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        Which_Fluid_In_Cell_U(nElemsX+ii,jj) = Which_Fluid_In_Cell_U(ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY    
      DO ii=1,nGhosts+1
        Which_Fluid_In_Cell_U(nElemsX+ii,jj) = Which_Fluid_In_Cell_U(nElemsX-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    fluid_status=INT( SIGN(1.0,PrimRefState2(nVar)) )
    IF (fluid_status .EQ. 0) THEN
      fluid_status=1
    END IF
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        Which_Fluid_In_Cell_U(nElemsX+ii,jj) = fluid_status
      END DO
    END DO
  CASE(4) ! Outflow
    fluid_status=INT( SIGN(1.0,PrimRefState2(nVar)) )
    IF (fluid_status .EQ. 0) THEN
      fluid_status=1
    END IF
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        Which_Fluid_In_Cell_U(nElemsX+ii,jj) = fluid_status
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=1,nGhosts+1
        Which_Fluid_In_Cell_U(nElemsX+ii,jj) = Which_Fluid_In_Cell_U(nElemsX-ii+1,jj)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition R not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!------------------------------!
! Top Boundary Conditions      !
!------------------------------!
SELECT CASE(BoundaryConditionsType(3))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        Which_Fluid_In_Cell_U(ii,nElemsY+jj) = Which_Fluid_In_Cell_U(ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        Which_Fluid_In_Cell_U(ii,nElemsY+jj) = Which_Fluid_In_Cell_U(ii,nElemsY-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    fluid_status=INT( SIGN(1.0,PrimRefState3(nVar)) )
    IF (fluid_status .EQ. 0) THEN
      fluid_status=1
    END IF
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        Which_Fluid_In_Cell_U(ii,nElemsY+jj) = fluid_status
      END DO
    END DO
  CASE(4) ! Outflow
    fluid_status=INT( SIGN(1.0,PrimRefState3(nVar)) )
    IF (fluid_status .EQ. 0) THEN
      fluid_status=1
    END IF
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        Which_Fluid_In_Cell_U(ii,nElemsY+jj) = fluid_status
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=1,nGhosts+1
        Which_Fluid_In_Cell_U(ii,nElemsY+jj) = Which_Fluid_In_Cell_U(ii,nElemsY-jj+1)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition T not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!------------------------------!
! Bottom Boundary Conditions   !
!------------------------------!
SELECT CASE(BoundaryConditionsType(1))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        Which_Fluid_In_Cell_U(ii,-nGhosts+jj) = Which_Fluid_In_Cell_U(ii,nElemsY-nGhosts+jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        Which_Fluid_In_Cell_U(ii,-nGhosts+jj) = Which_Fluid_In_Cell_U(ii,nGhosts-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    fluid_status=INT( SIGN(1.0,PrimRefState1(nVar)) )
    IF (fluid_status .EQ. 0) THEN
      fluid_status=1
    END IF
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        Which_Fluid_In_Cell_U(ii,-nGhosts+jj) = fluid_status
      END DO
    END DO
  CASE(4) ! Outflow
    fluid_status=INT( SIGN(1.0,PrimRefState1(nVar)) )
    IF (fluid_status .EQ. 0) THEN
      fluid_status=1
    END IF
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        Which_Fluid_In_Cell_U(ii,-nGhosts+jj) = fluid_status
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=0,nGhosts
        Which_Fluid_In_Cell_U(ii,-nGhosts+jj) = Which_Fluid_In_Cell_U(ii,nGhosts-jj+1)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition B not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!------------------------------!
! Left Corners Boundary Conditions!
!------------------------------!
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    !*Bottom left
    DO jj=-nGhosts,0
      DO ii=0,nGhosts
        Which_Fluid_In_Cell_U(-nGhosts+ii,jj) = Which_Fluid_In_Cell_U(nElemsX-nGhosts+ii,jj)
      END DO
    END DO

    !*Top left
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=0,nGhosts
        Which_Fluid_In_Cell_U(-nGhosts+ii,jj) = Which_Fluid_In_Cell_U(nElemsX-nGhosts+ii,jj)
      END DO
    END DO
END SELECT



!------------------------------!
! Right Corners Boundary Conditions    !
!------------------------------!
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    !*Bottom right
    DO jj=-nGhosts,0
      DO ii=1,nGhosts+1
        Which_Fluid_In_Cell_U(nElemsX+ii,jj) = Which_Fluid_In_Cell_U(ii,jj)
      END DO
    END DO

    !*Top right
    DO jj=nElemsY+1,nElemsY+1+nGhosts
      DO ii=1,nGhosts+1
        Which_Fluid_In_Cell_U(nElemsX+ii,jj) = Which_Fluid_In_Cell_U(ii,jj)
      END DO
    END DO
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE Impose_BC_on_Fluid_Cell_U
#endif
#endif
!===============================================================================!
!
!
!
!===============================================================================!
FUNCTION Switching_Function(x)
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: x
REAL            :: Switching_Function
! OLD
! REAL, PARAMETER :: x0=0.45
! REAL, PARAMETER :: alpha=2e-2

REAL, PARAMETER :: x0=0.15
REAL, PARAMETER :: x1=0.4
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!

!*I want 0 for eps=1
!*I want 1 for eps=0

! Switching_Function = EXP(-2000.0*x**6)
! Switching_Function = EXP(-3000.0*x**8)
! Switching_Function = EXP(-50000.0*x**10)

! IF (x .LE. x0) THEN
!   Switching_Function=0.0
! ELSE
!   Switching_Function=1.0-EXP(-((x - x0))**2/ alpha)
! END IF

  IF (x .LE. x0) THEN
     Switching_Function = 1.0-x**14 !*1.0-eps[mask1]**14
  ELSE IF ((x .GT. x0) .AND. (x .LT. x1)) THEN
     Switching_Function = EXP( 1.0 - 1.0 / (1.0 - ( ( (x - x0) / (x1-x0) ) ** 2) ) )*( (1.0-x0**14)-(1.0-x1)**14 )+(1.0-x1)**14
    !*result[mask2] =np.exp( 1.0 - 1.0 / (1.0 - ( ( (eps[mask2] - x0) / (x1-x0) ) ** 2) ) )*( (1.0-x0**14)-(1.0-x1)**14 )+(1.0-x1)**14
  ELSE
     Switching_Function = (1.0 - x)**14 !*(1-eps[mask3])**14
  END IF


! Switching_Function=0.0

!-------------------------------------------------------------------------------!
END FUNCTION Switching_Function
END MODULE MOD_Equation
!-------------------------------------------------------------------------------!
