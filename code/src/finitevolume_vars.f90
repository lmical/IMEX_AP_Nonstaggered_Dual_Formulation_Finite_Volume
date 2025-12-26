!===============================================================================!
MODULE MOD_FiniteVolume2D_vars
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PUBLIC
!-------------------------------------------------------------------------------!
! >> GLOBAL VARIABLES                                                           !
!-------------------------------------------------------------------------------!
#ifdef MULTIFLUID
INTEGER,PARAMETER   :: nVar  = 5
#else
INTEGER,PARAMETER   :: nVar  = 4
#endif
INTEGER,PARAMETER   :: nDims = 2
REAL                :: MESH_SX(1:nDims)
REAL                :: MESH_X0(1:nDims)
REAL                :: MESH_X1(1:nDims)
REAL                :: MESH_DX(1:nDims)
INTEGER             :: nGPs
INTEGER             :: nGhosts
INTEGER             :: nElemsX 
INTEGER             :: nElemsY

#ifdef PATANKAR
INTEGER             :: NNZsparse
INTEGER             :: NGlobalRows
INTEGER,ALLOCATABLE :: RowStart(:)
INTEGER,ALLOCATABLE :: ColumnsVector(:)
REAL, ALLOCATABLE   :: ProductionSparse(:)
REAL, ALLOCATABLE   :: DestructionSparse(:)
REAL,ALLOCATABLE    :: ProdUp(:,:)
REAL,ALLOCATABLE    :: DestUp(:,:)
INTEGER, ALLOCATABLE   :: SparseIndexMat(:,:)
#endif





REAL,ALLOCATABLE    :: MeshNodes(:,:,:)
REAL,ALLOCATABLE    :: MeshBary(:,:,:)
REAL,ALLOCATABLE    :: MeshGP(:,:,:,:,:)
REAL,ALLOCATABLE    :: WeightsGP(:,:)
REAL,ALLOCATABLE    :: WeightsGPBnd(:)
REAL                :: NormVectX(1:nDims)
REAL                :: NormVectY(1:nDims)
REAL                :: TangVectX(1:nDims)
REAL                :: TangVectY(1:nDims)


REAL,ALLOCATABLE    :: U(:,:,:)
REAL,ALLOCATABLE    :: V(:,:,:)
REAL,ALLOCATABLE    :: S(:,:,:)
REAL,ALLOCATABLE    :: SWB(:,:,:)
REAL,ALLOCATABLE    :: Ut(:,:,:)
REAL,ALLOCATABLE    :: FX(:,:,:)
REAL,ALLOCATABLE    :: FY(:,:,:)
REAL,ALLOCATABLE    :: FXWB(:,:,:)
REAL,ALLOCATABLE    :: FYWB(:,:,:)
REAL,ALLOCATABLE    :: WM(:,:,:,:)
REAL,ALLOCATABLE    :: WP(:,:,:,:)
REAL,ALLOCATABLE    :: FluxX(:,:,:,:)
REAL,ALLOCATABLE    :: FluxY(:,:,:,:)
LOGICAL,ALLOCATABLE :: Ind(:,:,:)

REAL,ALLOCATABLE    :: LLX(:,:,:,:)
REAL,ALLOCATABLE    :: RRX(:,:,:,:)
REAL,ALLOCATABLE    :: LLY(:,:,:,:)
REAL,ALLOCATABLE    :: RRY(:,:,:,:)
REAL,ALLOCATABLE    :: MCtoP(:,:,:,:)
REAL,ALLOCATABLE    :: MPtoC(:,:,:,:)


REAL,ALLOCATABLE    :: UN0(:,:,:)
REAL,ALLOCATABLE    :: K0(:,:,:)
REAL,ALLOCATABLE    :: K1(:,:,:)
REAL,ALLOCATABLE    :: K2(:,:,:)
REAL,ALLOCATABLE    :: K3(:,:,:)
REAL,ALLOCATABLE    :: K4(:,:,:)
REAL,ALLOCATABLE    :: K5(:,:,:)
REAL,ALLOCATABLE    :: FUp(:,:,:,:)
REAL,ALLOCATABLE    :: Up(:,:,:,:)
REAL,ALLOCATABLE    :: Ua(:,:,:,:)
INTEGER             :: MStepsMax

#ifdef CENTEREDPRIMITIVE
REAL,ALLOCATABLE    :: WC(:,:,:)
REAL,ALLOCATABLE    :: WCt(:,:,:)
REAL,ALLOCATABLE    :: UN0_WC(:,:,:)
REAL,ALLOCATABLE    :: K0_WC(:,:,:)
REAL,ALLOCATABLE    :: K1_WC(:,:,:)
REAL,ALLOCATABLE    :: K2_WC(:,:,:)
REAL,ALLOCATABLE    :: K3_WC(:,:,:)
REAL,ALLOCATABLE    :: K4_WC(:,:,:)
REAL,ALLOCATABLE    :: K5_WC(:,:,:)
REAL,ALLOCATABLE    :: FWp_WC(:,:,:,:)
REAL,ALLOCATABLE    :: Wp_WC(:,:,:,:)
REAL,ALLOCATABLE    :: Wa_WC(:,:,:,:)

!*The following structures are essential for IMEX, 
!*but I need them also for the explicit version, for the linear system to initialize the velocity in Gresho

INTEGER              :: NRows_WC
INTEGER, ALLOCATABLE :: From_2Indices_To_1Index_WC(:,:) !*ii,jj  -> ll     !*Also for ghosts
INTEGER, ALLOCATABLE :: From_1Index_To_2Indices_WC(:,:) !*ll,1:2 -> ii:jj  !*Also for ghosts
INTEGER, ALLOCATABLE :: Neighbours_1Index_WC(:,:)       !*ll,1:5 ->        !*Just for internal


#ifdef IMEX
!*-> For CRS of the matrix
INTEGER              :: NNZsparse_WC
REAL,    ALLOCATABLE :: Values_WC(:)
INTEGER, ALLOCATABLE :: Columns_WC(:)
INTEGER, ALLOCATABLE :: RowStart_WC(:)
REAL,    ALLOCATABLE :: Diagonal_WC(:)
!*-> For rhs of the system
REAL,    ALLOCATABLE :: rhs_WC(:)
!*-> For sol of the system
REAL,    ALLOCATABLE :: sol_WC(:)

INTEGER, ALLOCATABLE:: JacobiIterations_WC(:)
INTEGER             :: JacobiCounter_WC

#endif

#endif


#ifdef ACTIVEFLUX
!*For staggering in X direction
REAL,ALLOCATABLE    :: MeshBary_X(:,:,:)
REAL,ALLOCATABLE    :: MeshGP_X(:,:,:,:,:)
REAL,ALLOCATABLE    :: W_X(:,:,:)
REAL,ALLOCATABLE    :: Wt_X(:,:,:)
REAL,ALLOCATABLE    :: UN0_X(:,:,:)
REAL,ALLOCATABLE    :: K0_X(:,:,:)
REAL,ALLOCATABLE    :: K1_X(:,:,:)
REAL,ALLOCATABLE    :: K2_X(:,:,:)
REAL,ALLOCATABLE    :: K3_X(:,:,:)
REAL,ALLOCATABLE    :: K4_X(:,:,:)
REAL,ALLOCATABLE    :: K5_X(:,:,:)
REAL,ALLOCATABLE    :: FWp_X(:,:,:,:)
REAL,ALLOCATABLE    :: Wp_X(:,:,:,:)
REAL,ALLOCATABLE    :: Wa_X(:,:,:,:)

!*For staggering in Y direction
REAL,ALLOCATABLE    :: MeshBary_Y(:,:,:)
REAL,ALLOCATABLE    :: MeshGP_Y(:,:,:,:,:)
REAL,ALLOCATABLE    :: W_Y(:,:,:)
REAL,ALLOCATABLE    :: Wt_Y(:,:,:)
REAL,ALLOCATABLE    :: UN0_Y(:,:,:)
REAL,ALLOCATABLE    :: K0_Y(:,:,:)
REAL,ALLOCATABLE    :: K1_Y(:,:,:)
REAL,ALLOCATABLE    :: K2_Y(:,:,:)
REAL,ALLOCATABLE    :: K3_Y(:,:,:)
REAL,ALLOCATABLE    :: K4_Y(:,:,:)
REAL,ALLOCATABLE    :: K5_Y(:,:,:)
REAL,ALLOCATABLE    :: FWp_Y(:,:,:,:)
REAL,ALLOCATABLE    :: Wp_Y(:,:,:,:)
REAL,ALLOCATABLE    :: Wa_Y(:,:,:,:)

REAL,ALLOCATABLE    :: W_qp_X(:,:,:,:)
REAL,ALLOCATABLE    :: W_qp_Y(:,:,:,:)

!*The following structures are essential for IMEX, 
!*but I need them also for the explicit version, for the linear system to initialize the velocity in Gresho
!*In any case, they are referred to W_X and W_Y and so they are meant to be present if ACTIVEFLUX flag is active

!*Staggering in X direction
INTEGER              :: NRows_X
INTEGER, ALLOCATABLE :: From_2Indices_To_1Index_X(:,:) !*ii,jj  -> ll     !*Also for ghosts
INTEGER, ALLOCATABLE :: From_1Index_To_2Indices_X(:,:) !*ll,1:2 -> ii:jj  !*Also for ghosts
INTEGER, ALLOCATABLE :: Neighbours_1Index_X(:,:)       !*ll,1:5 ->        !*Just for internal

!*Staggering in Y direction
INTEGER              :: NRows_Y
INTEGER, ALLOCATABLE :: From_2Indices_To_1Index_Y(:,:) !*ii,jj  -> ll     !*Also for ghosts
INTEGER, ALLOCATABLE :: From_1Index_To_2Indices_Y(:,:) !*ll,1:2 -> ii:jj  !*Also for ghosts
INTEGER, ALLOCATABLE :: Neighbours_1Index_Y(:,:)       !*ll,1:5 ->        !*Just for internal

#if defined(IMEX) || defined(IMEXMOMENTUM)
!*For staggering in X direction
!*-> For CRS of the matrix
INTEGER              :: NNZsparse_X
REAL,    ALLOCATABLE :: Values_X(:)
INTEGER, ALLOCATABLE :: Columns_X(:)
INTEGER, ALLOCATABLE :: RowStart_X(:)
REAL,    ALLOCATABLE :: Diagonal_X(:)
!*-> For rhs of the system
REAL,    ALLOCATABLE :: rhs_X(:)
!*-> For sol of the system
REAL,    ALLOCATABLE :: sol_X(:)

INTEGER, ALLOCATABLE:: JacobiIterations_X(:)
INTEGER             :: JacobiCounter_X


!*For staggering in Y direction
INTEGER              :: NNZsparse_Y
REAL,    ALLOCATABLE :: Values_Y(:)
INTEGER, ALLOCATABLE :: Columns_Y(:)
INTEGER, ALLOCATABLE :: RowStart_Y(:)
REAL,    ALLOCATABLE :: Diagonal_Y(:)
!*-> For rhs of the system
REAL,    ALLOCATABLE :: rhs_Y(:)
!*-> For sol of the system
REAL,    ALLOCATABLE :: sol_Y(:)




INTEGER, ALLOCATABLE:: JacobiIterations_Y(:)
INTEGER             :: JacobiCounter_Y

#endif



#ifdef MULTIFLUID
INTEGER, ALLOCATABLE :: Troubled_Cell_U(:,:)
INTEGER, ALLOCATABLE :: Troubled_Cell_W_X(:,:)
INTEGER, ALLOCATABLE :: Troubled_Cell_W_Y(:,:)
INTEGER, ALLOCATABLE :: Which_Fluid_In_Cell_U(:,:)
#endif
#endif


#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
!*Shared
REAL,ALLOCATABLE    :: BX_Vol(:,:,:)
REAL,ALLOCATABLE    :: BY_Vol(:,:,:)
REAL,ALLOCATABLE    :: BX_Sur(:,:,:,:)     !*NB: Added a further dimension (with respect to FX and FY) corresponding to the specific side of the edge to take into account the wave speed
REAL,ALLOCATABLE    :: BY_Sur(:,:,:,:)     !*NB: Added a further dimension (with respect to FX and FY) corresponding to the specific side of the edge to take into account the wave speed
REAL,ALLOCATABLE    :: BX_Sur_qp(:,:,:,:)
REAL,ALLOCATABLE    :: BY_Sur_qp(:,:,:,:)
REAL,ALLOCATABLE    :: am_qp(:,:,:)        !*Local speed of propagation for each qp 
REAL,ALLOCATABLE    :: ap_qp(:,:,:)        !*Local speed of propagation for each qp 
#endif




REAL,ALLOCATABLE    :: UtWB(:,:,:)
REAL,ALLOCATABLE    :: Gravitational_Potential_Averages(:,:)

INTEGER,PARAMETER   :: UNIT_FILE = 123
INTEGER             :: WhichOutput
INTEGER             :: nOutputFiles
INTEGER             :: InitialCondition
INTEGER             :: BoundaryConditionsType(4)
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
INTEGER             :: BoundaryConditionsTypeIMEXLinearSystem(4)
#endif
#endif
REAL                :: PrimRefState1(1:nVar)
REAL                :: PrimRefState2(1:nVar)
REAL                :: PrimRefState3(1:nVar)
REAL                :: PrimRefState4(1:nVar)

REAL                :: t
REAL                :: tGlobal
REAL                :: dt
REAL                :: dt_Analyze
REAL                :: CFL
REAL                :: tEnd
REAL                :: Gmm
#ifdef MULTIFLUID
REAL                :: Gmm1
REAL                :: Gmm2
REAL                :: pinf1
REAL                :: pinf2
#endif
REAL                :: LambdaMaxX
REAL                :: LambdaMaxY
!*For very early time steps in particular situations
REAL                :: relaxed_dt
INTEGER             :: NRelaxedTimesteps
REAL                :: DiscontinuityLocation !*For Toro's Riemann Problems
REAL                :: Lambda_Lax
INTEGER             :: GravitationalPotentialFlag
#ifdef SW
REAL                :: Gravity
REAL                :: Kappa
#endif

INTEGER             :: maxTimeSteps
INTEGER             :: WhichRiemannSolver
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
INTEGER             :: WhichRiemannSolverPrimitiveSystem
#endif
REAL                :: computationalTime
REAL                :: global_min

#ifdef PATANKAR
INTEGER, ALLOCATABLE:: JacobiIterations(:)
INTEGER             :: JacobiCounter
#endif

#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#ifdef IMEX
REAL                :: max_ro, min_p
#endif
#endif

INTEGER             :: Reconstruction
INTEGER             :: ReconstructionFix
INTEGER             :: ReconstructedVariable
INTEGER             :: WhichSpeedEstimateForDtComputation
INTEGER             :: timescheme 
REAL,PARAMETER      :: wLobatto = 1./12.
REAL,PARAMETER      :: WENOEPS = 1.0E-6
INTEGER,PARAMETER   :: WENOEXP = 2 !*.0 !*ALERT
#ifdef ACTIVEFLUX
INTEGER             :: AF_PostProcessing
#endif

REAL                :: Ma_shock 
REAL                :: Ma_vortex
REAL                :: R_EOS


REAL,PARAMETER      :: PI           = ACOS(-1.0)
REAL,PARAMETER      :: EPS          = 1.0E-6
REAL,PARAMETER      :: ACCURACY     = 1.0E-30
REAL,PARAMETER      :: MIN_DENSITY  = 1.0E-6
REAL,PARAMETER      :: MIN_PRESSURE = 1.0E-30
REAL,PARAMETER      :: MIN_SPEED    = 0.0
REAL,PARAMETER      :: MIN_TIMESTEP = 1.0E-30
#ifdef RELAXATION
REAL                :: EPS_LM        !*LM=Low Mach or Lorenzo Micalizzi XD
#endif

REAL,PARAMETER      :: EPS_AVOIDING_DIVISION_BY_ZERO_IN_CU  = 1e-15 !*When the difference between the wave speeds is bigger than this parameter we switch to Rusanov
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#ifdef IMEX
REAL                :: SAFETY_K_TO_AVOID_HYPERBOLICITY_LOSS = 0.0!10.0
#endif
#ifdef MULTIFLUID
INTEGER, PARAMETER :: N_Troubled_Neighbors=1
#endif
#endif

CHARACTER(LEN=255)  :: VarNameVisu(1:nVar+1)
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
CHARACTER(LEN=255)  :: VarNameVisu_W(1:nVar+1)
#endif

#if defined(CENTEREDPRIMITIVE) && defined(PATHCONSERVATIVESHOCKDETECTION)
REAL                 :: K_coefficient_PCSD            = 0.025
INTEGER, PARAMETER   :: N_Limited_Neighbors_PCSD      = 0
INTEGER, ALLOCATABLE :: Cell_To_Limit_PCSD(:,:)
LOGICAL              :: IsSomeoneFlagged_PCSD         = .FALSE.
INTEGER              :: WhichVariableForLimiting_PCSD = 0 !*1 density, 2 momentum, 3 energy
#endif
#if defined(CENTEREDPRIMITIVE)
REAL,    ALLOCATABLE :: Errors_PCSD(:,:,:)
#endif


INTEGER             :: maxLenghtJacobiCounter

!*28
REAL, PARAMETER      :: Tau_Standard_Limiting   =0.5
REAL, PARAMETER      :: Theta_Standard_Limiting =2.0
!*29
REAL, PARAMETER      :: Tau_Overcompressive     =-0.25
REAL, PARAMETER      :: Theta_Overcompressive   =2.0

#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#ifdef IMEX
LOGICAL             :: HyperbolicityLoss = .FALSE.
#endif
#endif

!-------------------------------------------------------------------------------!
END MODULE MOD_FiniteVolume2D_vars
!===============================================================================!
