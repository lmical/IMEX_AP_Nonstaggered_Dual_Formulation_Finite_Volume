!===============================================================================!
MODULE MOD_FiniteVolume2D
!===============================================================================!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE InitializeFiniteVolume
  MODULE PROCEDURE InitializeFiniteVolume
END INTERFACE

INTERFACE FillInitialConditions
  MODULE PROCEDURE FillInitialConditions
END INTERFACE

INTERFACE InitializeWBVariables
  MODULE PROCEDURE InitializeWBVariables
END INTERFACE

INTERFACE FVTimeDerivative
  MODULE PROCEDURE FVTimeDerivative
END INTERFACE

INTERFACE FinalizeFiniteVolume
  MODULE PROCEDURE FinalizeFiniteVolume
END INTERFACE

INTERFACE ComputeTransitionMatricesConservedInput
  MODULE PROCEDURE ComputeTransitionMatricesConservedInput
END INTERFACE

INTERFACE Impose_Symmetric_Update_Y_Axis
  MODULE PROCEDURE Impose_Symmetric_Update_Y_Axis
END INTERFACE

#ifdef CENTEREDPRIMITIVE
INTERFACE Overwrite_WC_from_U
  MODULE PROCEDURE Overwrite_WC_from_U
END INTERFACE

#ifdef PATHCONSERVATIVESHOCKDETECTION
INTERFACE Shock_Detector_Based_On_Path_Conservative
  MODULE PROCEDURE Shock_Detector_Based_On_Path_Conservative
END INTERFACE
#endif
#endif



#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
INTERFACE FVTimeDerivativePrimitiveOnly
  MODULE PROCEDURE FVTimeDerivativePrimitiveOnly
END INTERFACE

INTERFACE FVTimeDerivativeConservedOnly
  MODULE PROCEDURE FVTimeDerivativeConservedOnly
END INTERFACE
#endif

#if defined(CENTEREDPRIMITIVE)
INTERFACE Compute_Errors_PCSD
  MODULE PROCEDURE Compute_Errors_PCSD
END INTERFACE
#endif

#ifdef ACTIVEFLUX
INTERFACE AF_PostProcessing_Subroutine
  MODULE PROCEDURE AF_PostProcessing_Subroutine
END INTERFACE

#if defined(IMEX) || defined(IMEXMOMENTUM)
INTERFACE Put_Matrix_In_Vector_X
  MODULE PROCEDURE Put_Matrix_In_Vector_X
END INTERFACE

INTERFACE Put_Vector_In_Matrix_X
  MODULE PROCEDURE Put_Vector_In_Matrix_X
END INTERFACE

INTERFACE Put_Matrix_In_Vector_Y
  MODULE PROCEDURE Put_Matrix_In_Vector_Y
END INTERFACE

INTERFACE Put_Vector_In_Matrix_Y
  MODULE PROCEDURE Put_Vector_In_Matrix_Y
END INTERFACE
#endif
#endif

#ifdef CENTEREDPRIMITIVE
#if defined(IMEX) || defined(IMEXMOMENTUM)
INTERFACE Put_Matrix_In_Vector
  MODULE PROCEDURE Put_Matrix_In_Vector
END INTERFACE

INTERFACE Put_Vector_In_Matrix
  MODULE PROCEDURE Put_Vector_In_Matrix
END INTERFACE
#endif
#endif

!-------------------------------------------------------------------------------!
PUBLIC :: InitializeFiniteVolume
PUBLIC :: FillInitialConditions
PUBLIC :: InitializeWBVariables
PUBLIC :: FVTimeDerivative
PUBLIC :: FinalizeFiniteVolume
PUBLIC :: ComputeTransitionMatricesConservedInput
PUBLIC :: Impose_Symmetric_Update_Y_Axis
#ifdef CENTEREDPRIMITIVE
PUBLIC :: Overwrite_WC_from_U
#ifdef PATHCONSERVATIVESHOCKDETECTION
PUBLIC :: Shock_Detector_Based_On_Path_Conservative
#endif
#endif
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
PUBLIC :: FVTimeDerivativePrimitiveOnly
PUBLIC :: FVTimeDerivativeConservedOnly
#endif
#if defined(CENTEREDPRIMITIVE)
PUBLIC :: Compute_Errors_PCSD
#endif
#ifdef ACTIVEFLUX
PUBLIC :: AF_PostProcessing_Subroutine
#if defined(IMEX) || defined(IMEXMOMENTUM)
PUBLIC :: Put_Matrix_In_Vector_X
PUBLIC :: Put_Vector_In_Matrix_X
PUBLIC :: Put_Matrix_In_Vector_Y
PUBLIC :: Put_Vector_In_Matrix_Y
#endif
#endif
#ifdef CENTEREDPRIMITIVE
#if defined(IMEX) || defined(IMEXMOMENTUM)
PUBLIC :: Put_Matrix_In_Vector
PUBLIC :: Put_Vector_In_Matrix
#endif
#endif
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
SUBROUTINE InitializeFiniteVolume()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MeshNodes
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP  
USE MOD_FiniteVolume2D_vars,ONLY: NormVectX
USE MOD_FiniteVolume2D_vars,ONLY: NormVectY
USE MOD_FiniteVolume2D_vars,ONLY: TangVectX
USE MOD_FiniteVolume2D_vars,ONLY: TangVectY
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: UtWB
USE MOD_FiniteVolume2D_vars,ONLY: S
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: FX
USE MOD_FiniteVolume2D_vars,ONLY: FY
USE MOD_FiniteVolume2D_vars,ONLY: FXWB
USE MOD_FiniteVolume2D_vars,ONLY: FYWB
USE MOD_FiniteVolume2D_vars,ONLY: SWB
USE MOD_FiniteVolume2D_vars,ONLY: FluxX
USE MOD_FiniteVolume2D_vars,ONLY: FluxY
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction

USE MOD_FiniteVolume2D_vars,ONLY: UN0
USE MOD_FiniteVolume2D_vars,ONLY: K0
USE MOD_FiniteVolume2D_vars,ONLY: K1
USE MOD_FiniteVolume2D_vars,ONLY: K2
USE MOD_FiniteVolume2D_vars,ONLY: K3
USE MOD_FiniteVolume2D_vars,ONLY: K4
USE MOD_FiniteVolume2D_vars,ONLY: K5
USE MOD_FiniteVolume2D_vars,ONLY: Up
USE MOD_FiniteVolume2D_vars,ONLY: Ua
USE MOD_FiniteVolume2D_vars,ONLY: FUp
USE MOD_FiniteVolume2D_vars,ONLY: Gravitational_Potential_Averages
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGP
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGPBnd
USE MOD_FiniteVolume2D_vars,ONLY: timescheme 
USE MOD_FiniteVolume2D_vars,ONLY: maxTimeSteps 
USE MOD_FiniteVolume2D_vars,ONLY: MStepsMax


USE MOD_FiniteVolume2D_vars,ONLY: LLX
USE MOD_FiniteVolume2D_vars,ONLY: RRX
USE MOD_FiniteVolume2D_vars,ONLY: LLY
USE MOD_FiniteVolume2D_vars,ONLY: RRY
USE MOD_FiniteVolume2D_vars,ONLY: MCtoP
USE MOD_FiniteVolume2D_vars,ONLY: MPtoC


#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: UN0_WC
USE MOD_FiniteVolume2D_vars,ONLY: K0_WC
USE MOD_FiniteVolume2D_vars,ONLY: K1_WC
USE MOD_FiniteVolume2D_vars,ONLY: K2_WC
USE MOD_FiniteVolume2D_vars,ONLY: K3_WC
USE MOD_FiniteVolume2D_vars,ONLY: K4_WC
USE MOD_FiniteVolume2D_vars,ONLY: K5_WC
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC


!*The following structures are essential for IMEX, 
!*but I need them also for the explicit version, for the linear system to initialize the velocity in Gresho

USE MOD_FiniteVolume2D_vars,ONLY: NRows_WC
USE MOD_FiniteVolume2D_vars,ONLY: From_2Indices_To_1Index_WC
USE MOD_FiniteVolume2D_vars,ONLY: From_1Index_To_2Indices_WC
USE MOD_FiniteVolume2D_vars,ONLY: Neighbours_1Index_WC

#ifdef IMEX
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_WC
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_WC
USE MOD_FiniteVolume2D_vars,ONLY: Values_WC
USE MOD_FiniteVolume2D_vars,ONLY: Columns_WC
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_WC
USE MOD_FiniteVolume2D_vars,ONLY: rhs_WC
USE MOD_FiniteVolume2D_vars,ONLY: sol_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_WC
#endif
#endif

#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
#if defined(IMEX) || defined(IMEXMOMENTUM)
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsTypeIMEXLinearSystem
USE MOD_FiniteVolume2D_vars,ONLY: maxLenghtJacobiCounter
#endif
#endif


#if defined(CENTEREDPRIMITIVE) && defined(PATHCONSERVATIVESHOCKDETECTION)
USE MOD_FiniteVolume2D_vars,ONLY: Cell_To_Limit_PCSD
#endif
#if defined(CENTEREDPRIMITIVE)
USE MOD_FiniteVolume2D_vars,ONLY: Errors_PCSD
#endif



#ifdef ACTIVEFLUX
!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary_X
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_X  
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X
USE MOD_FiniteVolume2D_vars,ONLY: UN0_X
USE MOD_FiniteVolume2D_vars,ONLY: K0_X
USE MOD_FiniteVolume2D_vars,ONLY: K1_X
USE MOD_FiniteVolume2D_vars,ONLY: K2_X
USE MOD_FiniteVolume2D_vars,ONLY: K3_X
USE MOD_FiniteVolume2D_vars,ONLY: K4_X
USE MOD_FiniteVolume2D_vars,ONLY: K5_X
USE MOD_FiniteVolume2D_vars,ONLY: Wp_X
USE MOD_FiniteVolume2D_vars,ONLY: Wa_X
USE MOD_FiniteVolume2D_vars,ONLY: FWp_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary_Y
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_Y
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y
USE MOD_FiniteVolume2D_vars,ONLY: UN0_Y
USE MOD_FiniteVolume2D_vars,ONLY: K0_Y
USE MOD_FiniteVolume2D_vars,ONLY: K1_Y
USE MOD_FiniteVolume2D_vars,ONLY: K2_Y
USE MOD_FiniteVolume2D_vars,ONLY: K3_Y
USE MOD_FiniteVolume2D_vars,ONLY: K4_Y
USE MOD_FiniteVolume2D_vars,ONLY: K5_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wp_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wa_Y
USE MOD_FiniteVolume2D_vars,ONLY: FWp_Y

USE MOD_FiniteVolume2D_vars,ONLY: W_qp_X
USE MOD_FiniteVolume2D_vars,ONLY: W_qp_Y

!*The following structures are essential for IMEX, 
!*but I need them also for the explicit version, for the linear system to initialize the velocity in Gresho
!*In any case, they are referred to W_X and W_Y and so they are meant to be present if ACTIVEFLUX flag is active

!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: NRows_X
USE MOD_FiniteVolume2D_vars,ONLY: From_2Indices_To_1Index_X
USE MOD_FiniteVolume2D_vars,ONLY: From_1Index_To_2Indices_X
USE MOD_FiniteVolume2D_vars,ONLY: Neighbours_1Index_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: NRows_Y
USE MOD_FiniteVolume2D_vars,ONLY: From_2Indices_To_1Index_Y
USE MOD_FiniteVolume2D_vars,ONLY: From_1Index_To_2Indices_Y
USE MOD_FiniteVolume2D_vars,ONLY: Neighbours_1Index_Y

#if defined(IMEX) || defined(IMEXMOMENTUM)
!*For staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_X
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_X
USE MOD_FiniteVolume2D_vars,ONLY: Values_X
USE MOD_FiniteVolume2D_vars,ONLY: Columns_X
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_X
USE MOD_FiniteVolume2D_vars,ONLY: rhs_X
USE MOD_FiniteVolume2D_vars,ONLY: sol_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_X

!*For staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse_Y
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_Y
USE MOD_FiniteVolume2D_vars,ONLY: Values_Y
USE MOD_FiniteVolume2D_vars,ONLY: Columns_Y
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_Y
USE MOD_FiniteVolume2D_vars,ONLY: rhs_Y
USE MOD_FiniteVolume2D_vars,ONLY: sol_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_Y
#endif

#ifdef MULTIFLUID
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_U
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_W_X
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Which_Fluid_In_Cell_U
#endif

#endif

#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
!*Shared
USE MOD_FiniteVolume2D_vars,ONLY: BX_Vol
USE MOD_FiniteVolume2D_vars,ONLY: BY_Vol
USE MOD_FiniteVolume2D_vars,ONLY: BX_Sur
USE MOD_FiniteVolume2D_vars,ONLY: BY_Sur
USE MOD_FiniteVolume2D_vars,ONLY: BX_Sur_qp
USE MOD_FiniteVolume2D_vars,ONLY: BY_Sur_qp
USE MOD_FiniteVolume2D_vars,ONLY: am_qp
USE MOD_FiniteVolume2D_vars,ONLY: ap_qp
#endif



#ifdef PATANKAR 
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse
USE MOD_FiniteVolume2D_vars,ONLY: NGlobalRows
USE MOD_FiniteVolume2D_vars,ONLY: ProductionSparse
USE MOD_FiniteVolume2D_vars,ONLY: DestructionSparse
USE MOD_FiniteVolume2D_vars,ONLY: RowStart
USE MOD_FiniteVolume2D_vars,ONLY: ColumnsVector
USE MOD_FiniteVolume2D_vars,ONLY: ProdUp
USE MOD_FiniteVolume2D_vars,ONLY: DestUp 
USE MOD_FiniteVolume2D_vars,ONLY: SparseIndexMat
USE MOD_Mesh,               ONLY: GlobalElem
USE MOD_FiniteVolume2D_vars,ONLY: JacobiCounter
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, iSparseVect, iRowStart, K, L
CHARACTER(LEN=255) :: ErrorMessage
INTEGER :: indc
!-------------------------------------------------------------------------------!

SELECT CASE (Reconstruction)
#ifndef ACTIVEFLUX
  !*IF NOT ACTIVE FLUX RECONSTRUCTION IS POSITIVE
  CASE(1,10) ! NONE
    nGhosts = 1
    nGPs    = 1
  CASE(2,20,21,22,23,24,25,26,27,28,29) ! MUSCL
    nGhosts = 1
    nGPs    = 1
  CASE(3) ! WENO3
    nGhosts = 1
    nGPs    = 2
  CASE(4,5) ! WENO5
    nGhosts = 2
    nGPs    = 4  !*NB: In principle 3 nGPs would have been ok but there were negative weights
  CASE(7) ! WENO7
    nGhosts = 3
    nGPs    = 4
#else
  !*IF ACTIVE FLUX RECONSTRUCTION IS NEGATIVE
  CASE(-1,-10) ! NONE
    nGhosts = 1
    nGPs    = 1
  CASE(-2,-20,-21,-22,-23,-24,-25,-26,-27,-28,-29) ! MUSCL
    nGhosts = 1
    nGPs    = 1
  CASE(-3) ! WENO3
    nGhosts = 1
    nGPs    = 2
  CASE(-4,-5) ! WENO5
    nGhosts = 2
    nGPs    = 4  !*NB: In principle 3 nGPs would have been ok but there were negative weights
  CASE(-7) ! WENO7
    nGhosts = 3
    nGPs    = 4
#endif
  CASE DEFAULT
    ErrorMessage = "Reconstruction not implemented. NB: Check the sign. It may be my trick to distinguish ACTIVE FLUX and NOT ACTIVE FLUX results"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

ALLOCATE(WeightsGP(1:nGPs, 1:nGPs))
ALLOCATE(WeightsGPBnd(1:nGPs))

ALLOCATE(MeshNodes(1:nDims,       0:nElemsX,0:nElemsY))
ALLOCATE(MeshBary (1:nDims,       1:nElemsX,1:nElemsY))
ALLOCATE(MeshGP   (1:nDims,       1:nElemsX,1:nElemsY,1:nGPs,1:nGPs))



ALLOCATE( U(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
ALLOCATE( V(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
ALLOCATE(Ut(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(Gravitational_Potential_Averages(1:nElemsX,1:nElemsY))
#ifdef MULTIFLUID
!*Adding one component for PWC reconstruction of gmm and pinf
ALLOCATE(WM(1:nVar+1,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)) !*NB: Adding 1 in both dimensions for staggered meshes (used only in ACTIVE FLUX)
ALLOCATE(WP(1:nVar+1,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)) !*NB: Adding 1 in both dimensions for staggered meshes (used only in ACTIVE FLUX)
#else
ALLOCATE(WM(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)) 
ALLOCATE(WP(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)) 
#endif
ALLOCATE( S(1:nVar,1:nElemsX+1,1:nElemsY+1))            !*NB: Adding 1 in both dimensions for staggered meshes (used only in ACTIVE FLUX)
ALLOCATE(FX(1:nVar,0:nElemsX+1,1:nElemsY+1))            !*NB: Adding 1 in both dimensions for staggered meshes (used only in ACTIVE FLUX)
ALLOCATE(FY(1:nVar,1:nElemsX+1,0:nElemsY+1))            !*NB: Adding 1 in both dimensions for staggered meshes (used only in ACTIVE FLUX)
ALLOCATE(FluxX(1:nVar,1:nGPs,0:nElemsX+1,1:nElemsY+1))  !*NB: Adding 1 in both dimensions for staggered meshes (used only in ACTIVE FLUX)
ALLOCATE(FluxY(1:nVar,1:nGPs,1:nElemsX+1,0:nElemsY+1))  !*NB: Adding 1 in both dimensions for staggered meshes (used only in ACTIVE FLUX)
ALLOCATE(Ind(1:2,0:nElemsX+1,0:nElemsY+1))




ALLOCATE(UN0(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K0(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K1(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K2(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K3(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K4(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K5(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(Ua(1:MStepsMax,1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(Up(1:MStepsMax,1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(FUp(1:MStepsMax,1:nVar,1:nElemsX,1:nElemsY))

!*NB: Allocated with +1 in X and Y direction
ALLOCATE(LLX( 1:nVar,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1+1) ) 
ALLOCATE(RRX( 1:nVar,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1+1) ) 
ALLOCATE(LLY( 1:nVar,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1+1) ) 
ALLOCATE(RRY( 1:nVar,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1+1) ) 
ALLOCATE(MCtoP( 1:nVar,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1+1) ) 
ALLOCATE(MPtoC( 1:nVar,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1+1) ) 



#ifdef CENTEREDPRIMITIVE
ALLOCATE( WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)) !*Like standard
ALLOCATE( WCt(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)) !*NB:made analogous to WC for IMEX
ALLOCATE(UN0_WC(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K0_WC(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K1_WC(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K2_WC(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K3_WC(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K4_WC(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(K5_WC(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(Wa_WC(1:MStepsMax,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))!*NB:made analogous to WC for IMEX
ALLOCATE(Wp_WC(1:MStepsMax,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))!*NB:made analogous to WC for IMEX
ALLOCATE(FWp_WC(1:MStepsMax,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))!*NB:made analogous to WC for IMEX
#endif

#if defined(CENTEREDPRIMITIVE) && defined(PATHCONSERVATIVESHOCKDETECTION)
ALLOCATE( Cell_To_Limit_PCSD(-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)) 
#endif
#if defined(CENTEREDPRIMITIVE)
ALLOCATE( Errors_PCSD(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)) 
#endif


#ifdef ACTIVEFLUX
!*Staggering in X direction
!*NB: +1 for X
ALLOCATE(MeshBary_X (1:nDims,       1:nElemsX+1,1:nElemsY))                   !*NB:+1 in X direction with respect to MeshBary
ALLOCATE(MeshGP_X   (1:nDims,       1:nElemsX+1,1:nElemsY,1:nGPs,1:nGPs))     !*NB:+1 in X direction with respect to MeshGP  
ALLOCATE(W_X( 1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1))!*NB:+1 in X direction with respect to U
ALLOCATE(Wt_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1))!*NB:made analogous to W_X for IMEX

ALLOCATE(UN0_X(1:nVar,1:nElemsX+1,1:nElemsY))                                 !*NB:+1 in X direction with respect to original counterpart
ALLOCATE(K0_X(1:nVar,1:nElemsX+1,1:nElemsY))                                  !*NB:+1 in X direction with respect to original counterpart 
ALLOCATE(K1_X(1:nVar,1:nElemsX+1,1:nElemsY))                                  !*NB:+1 in X direction with respect to original counterpart 
ALLOCATE(K2_X(1:nVar,1:nElemsX+1,1:nElemsY))                                  !*NB:+1 in X direction with respect to original counterpart 
ALLOCATE(K3_X(1:nVar,1:nElemsX+1,1:nElemsY))                                  !*NB:+1 in X direction with respect to original counterpart   
ALLOCATE(K4_X(1:nVar,1:nElemsX+1,1:nElemsY))                                  !*NB:+1 in X direction with respect to original counterpart 
ALLOCATE(K5_X(1:nVar,1:nElemsX+1,1:nElemsY))                                  !*NB:+1 in X direction with respect to original counterpart 
ALLOCATE(Wa_X( 1:MStepsMax,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1))!*NB:made analogous to W_X for IMEX
ALLOCATE(Wp_X( 1:MStepsMax,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1))!*NB:made analogous to W_X for IMEX
ALLOCATE(FWp_X(1:MStepsMax,1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1))!*NB:made analogous to W_X for IMEX


!*Staggering in Y direction
!*NB: +1 for Y
ALLOCATE(MeshBary_Y (1:nDims,       1:nElemsX,1:nElemsY+1))                   !*NB:+1 in Y direction with respect to MeshBary
ALLOCATE(MeshGP_Y   (1:nDims,       1:nElemsX,1:nElemsY+1,1:nGPs,1:nGPs))     !*NB:+1 in Y direction with respect to MeshGP  
ALLOCATE(W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)) !*NB:+1 in Y direction with respect to U
ALLOCATE(Wt_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1))!*NB:made analogous to W_Y for IMEX

ALLOCATE(UN0_Y(1:nVar,1:nElemsX,1:nElemsY+1))                                 !*NB:+1 in Y direction with respect to original counterpart
ALLOCATE(K0_Y(1:nVar,1:nElemsX,1:nElemsY+1))                                  !*NB:+1 in Y direction with respect to original counterpart 
ALLOCATE(K1_Y(1:nVar,1:nElemsX,1:nElemsY+1))                                  !*NB:+1 in Y direction with respect to original counterpart 
ALLOCATE(K2_Y(1:nVar,1:nElemsX,1:nElemsY+1))                                  !*NB:+1 in Y direction with respect to original counterpart 
ALLOCATE(K3_Y(1:nVar,1:nElemsX,1:nElemsY+1))                                  !*NB:+1 in Y direction with respect to original counterpart   
ALLOCATE(K4_Y(1:nVar,1:nElemsX,1:nElemsY+1))                                  !*NB:+1 in Y direction with respect to original counterpart 
ALLOCATE(K5_Y(1:nVar,1:nElemsX,1:nElemsY+1))                                  !*NB:+1 in Y direction with respect to original counterpart 
ALLOCATE(Wa_Y( 1:MStepsMax,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1))!*NB:made analogous to W_Y for IMEX
ALLOCATE(Wp_Y( 1:MStepsMax,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1))!*NB:made analogous to W_Y for IMEX
ALLOCATE(FWp_Y(1:MStepsMax,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1))!*NB:made analogous to W_Y for IMEX


ALLOCATE(W_qp_X(1:nVar,1:nGPs,0:nElemsX,1:nElemsY))                           !*NB: Same dimensions as (original) FluxX 
ALLOCATE(W_qp_Y(1:nVar,1:nGPs,1:nElemsX,0:nElemsY))                           !*NB: Same dimensions as (original) FluxY 


#ifdef MULTIFLUID
ALLOCATE( Troubled_Cell_U(-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)) 
ALLOCATE(Troubled_Cell_W_X(-nGhosts:nElemsX+nGhosts+1+1+1,-nGhosts:nElemsY+nGhosts+1))!*NB:+1 in X direction with respect to U
ALLOCATE(Troubled_Cell_W_Y(-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1+1))!*NB:+1 in Y direction with respect to U
ALLOCATE( Which_Fluid_In_Cell_U(-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)) 
#endif


#endif

#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
!*Shared
ALLOCATE(BX_Vol(1:nVar,1:nElemsX+1,1:nElemsY+1))                              !*NB: Like S but both dimensions+1
ALLOCATE(BY_Vol(1:nVar,1:nElemsX+1,1:nElemsY+1))                              !*NB: Like S but both dimensions+1
ALLOCATE(BX_Sur(1:nVar,0:nElemsX+1,1:nElemsY+1,1:2))                          !*NB: Like FX    (with +1 in both directions) !*NB: Added a further dimension corresponding to the specific side of the edge to take into account the wave speed
ALLOCATE(BY_Sur(1:nVar,1:nElemsX+1,0:nElemsY+1,1:2))                          !*NB: Like FY    (with +1 in both directions) !*NB: Added a further dimension corresponding to the specific side of the edge to take into account the wave speed
ALLOCATE(BX_Sur_qp(1:nVar,1:nGPs,0:nElemsX+1,1:nElemsY+1))                    !*NB: Like FluxX (with +1 in both directions)
ALLOCATE(BY_Sur_qp(1:nVar,1:nGPs,1:nElemsX+1,0:nElemsY+1))                    !*NB: Like FluxX (with +1 in both directions)
ALLOCATE(am_qp(1:nGPs,0:nElemsX+1,0:nElemsY+1))                               !*NB: To be used for BX_Sur_qp and BY_Sur_qp, so the allocation starts from 0 in both X and Y dimensions
ALLOCATE(ap_qp(1:nGPs,0:nElemsX+1,0:nElemsY+1))                               !*NB: To be used for BX_Sur_qp and BY_Sur_qp, so the allocation starts from 0 in both X and Y dimensions
#endif


#ifdef WELLBALANCED
ALLOCATE( UtWB(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE( SWB(1:nVar,1:nElemsX,1:nElemsY))
ALLOCATE(FXWB(1:nVar,0:nElemsX,1:nElemsY))
ALLOCATE(FYWB(1:nVar,1:nElemsX,0:nElemsY))
#endif

#ifdef PATANKAR
NNZsparse = 5*nElemsX*nElemsY
NGlobalRows = nElemsX*nElemsY
ALLOCATE(ProductionSparse(1:NNZsparse))
ALLOCATE(DestructionSparse(1:NNZsparse))
ALLOCATE(ColumnsVector(1:NNZsparse))
ALLOCATE(RowStart(1:NGlobalRows+1))
ALLOCATE(ProdUp(1:4,1:NNZsparse))
ALLOCATE(DestUp(1:4,1:NNZsparse))
ALLOCATE(SparseIndexMat(1:nElemsX*nElemsY,1:5))
ALLOCATE(JacobiIterations(1:maxTimeSteps))
#endif

U  = 0.0
V  = 0.0
S  = 0.0
Ut = 0.0
FX = 0.0
FY = 0.0
WM = 0.0
WP = 0.0
FluxX = 0.0
FluxY = 0.0


Ind = .FALSE.

UN0 = 0.0
K0  = 0.0
K1  = 0.0
K2  = 0.0
K3  = 0.0
K4  = 0.0
K5  = 0.0

Gravitational_Potential_Averages = 0.

Ua  = 0.
Up  = 0.
FUp = 0.

#ifdef CENTEREDPRIMITIVE
WC     = 0.0
WCt    = 0.0
UN0_WC = 0.0
K0_WC  = 0.0
K1_WC  = 0.0
K2_WC  = 0.0
K3_WC  = 0.0
K4_WC  = 0.0
K5_WC  = 0.0
Wa_WC  = 0.0
Wp_WC  = 0.0
FWp_WC = 0.0
#endif

#if defined(CENTEREDPRIMITIVE) && defined(PATHCONSERVATIVESHOCKDETECTION)
Cell_To_Limit_PCSD=0 
Errors_PCSD       =0.0
#endif


#ifdef ACTIVEFLUX
!*Staggering in X direction
W_X=0.
UN0_X = 0.0
K0_X  = 0.0
K1_X  = 0.0
K2_X  = 0.0
K3_X  = 0.0
K4_X  = 0.0
K5_X  = 0.0
Wa_X  = 0.0
Wp_X  = 0.0
FWp_X = 0.0

!*Staggering in Y direction
W_Y=0.
UN0_Y = 0.0
K0_Y  = 0.0
K1_Y  = 0.0
K2_Y  = 0.0
K3_Y  = 0.0
K4_Y  = 0.0
K5_Y  = 0.0
Wa_Y  = 0.0
Wp_Y  = 0.0
FWp_Y = 0.0

!*Shared
BX_Vol    = 0.
BY_Vol    = 0.
BX_Sur    = 0.
BY_Sur    = 0.
BX_Sur_qp = 0.
BY_Sur_qp = 0.
am_qp     = 0.
ap_qp     = 0.

W_qp_X = 0.
W_qp_Y = 0.


#ifdef MULTIFLUID
Troubled_Cell_U=0
Troubled_Cell_W_X=0
Troubled_Cell_W_Y=0
#endif

#endif




#ifdef WELLBALANCED
UtWB = 0.
SWB  = 0.0
FXWB = 0.0
FYWB = 0.0
#endif



#ifdef PATANKAR
iSparseVect = 0
iRowStart=0

DO jj=1,nElemsY
  DO ii=1,nElemsX
    K = GlobalElem(ii,jj)
    !K = (ii-1)*NelemsY + jj
    iSparseVect = iSparseVect + 1
    iRowStart = iRowStart +1

    RowStart(iRowStart) = iSparseVect
    !*PRINT*, "iRowStart", iRowStart, "RowStart(iRowStart)", RowStart(iRowStart)


    L = GlobalElem(ii,jj-1)
    ColumnsVector(iSparseVect) = L
    SparseIndexMat(K,1) = iSparseVect


    iSparseVect = iSparseVect + 1
    L = GlobalElem(ii-1,jj)
    ColumnsVector(iSparseVect) = L
    SparseIndexMat(K,2) = iSparseVect


    iSparseVect = iSparseVect + 1
    L = GlobalElem(ii,jj)
    ColumnsVector(iSparseVect) = L
    SparseIndexMat(K,3) = iSparseVect


    iSparseVect = iSparseVect + 1
    L = GlobalElem(ii+1,jj)
    ColumnsVector(iSparseVect) = L
    SparseIndexMat(K,4) = iSparseVect



    iSparseVect = iSparseVect + 1
    L = GlobalElem(ii,jj+1)
    ColumnsVector(iSparseVect) = L
    SparseIndexMat(K,5) = iSparseVect


  END DO
END DO

RowStart(nElemsX*nElemsY+1) = RowStart(nElemsX*nElemsY)+5 
!*PRINT*, "iRowStart", nElemsX*nElemsY+1, "RowStart(iRowStart)", RowStart(nElemsX*nElemsY+1)

JacobiCounter = 0

#endif


!*=================================================================
!* IMEX FOR CENTEREDPRIMITIVE
!*=================================================================
!*Structures for IMEX linear system
!*NB: The following is kind of analogous to what previously done
!*But it is organized slightly differently and more general

!*We only care about the first layer
!*The other one will have -1 (non-physical) values in our structures


#ifdef CENTEREDPRIMITIVE
!*The following structures are essential for IMEX, 
!*but I need them also for the explicit version, for the linear system to initialize the velocity in Gresho


NRows_WC=nElemsX*nElemsY+2*nElemsX+2*nElemsY
ALLOCATE(From_2Indices_To_1Index_WC(-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)) !*ii,jj  -> ll     !*Also for ghosts, same as U !*But careful with corners (-1)
ALLOCATE(From_1Index_To_2Indices_WC(1:NRows_WC,1:2))                                           !*ll,1:2 -> ii:jj  !*Also for ghosts excluding corners
ALLOCATE(Neighbours_1Index_WC(1:(nElemsX*nElemsY),1:5))                                       !*ll,1:5 ->        !*Just for internal

From_2Indices_To_1Index_WC = -1 !*Initialization to -1 meaning non-physical
From_1Index_To_2Indices_WC = -1 !*Initialization to -1 meaning non-physical
Neighbours_1Index_WC       = -1 !*Initialization to -1 meaning non-physical

!*-------------------------
!*Numeration
!*-------------------------

indc=0
!*->Inside the domain
DO jj=1,nElemsY
  DO ii=1,nElemsX 
    indc=indc+1
    From_2Indices_To_1Index_WC(ii,jj)=indc
    ! PRINT*, ii,jj,"->",indc
    From_1Index_To_2Indices_WC(indc,1)=ii
    From_1Index_To_2Indices_WC(indc,2)=jj

  END DO
END DO

!*->Down boundary
jj=0
DO ii=1,nElemsX
  indc=indc+1
  From_2Indices_To_1Index_WC(ii,jj)=indc
  ! PRINT*, ii,jj,"->",indc
  From_1Index_To_2Indices_WC(indc,1)=ii
  From_1Index_To_2Indices_WC(indc,2)=jj
END DO

!*->Right boundary
ii=nElemsX+1
DO jj=1,nElemsY
  indc=indc+1
  From_2Indices_To_1Index_WC(ii,jj)=indc
  ! PRINT*, ii,jj,"->",indc
  From_1Index_To_2Indices_WC(indc,1)=ii
  From_1Index_To_2Indices_WC(indc,2)=jj
END DO

!*->Top boundary
jj=nElemsY+1
DO ii=1,nElemsX 
  indc=indc+1
  From_2Indices_To_1Index_WC(ii,jj)=indc
  ! PRINT*, ii,jj,"->",indc
  From_1Index_To_2Indices_WC(indc,1)=ii
  From_1Index_To_2Indices_WC(indc,2)=jj
END DO

!*->Left boundary
ii=0
DO jj=1,nElemsY
  indc=indc+1
  From_2Indices_To_1Index_WC(ii,jj)=indc
  ! PRINT*, ii,jj,"->",indc
  From_1Index_To_2Indices_WC(indc,1)=ii
  From_1Index_To_2Indices_WC(indc,2)=jj
END DO

!*-----------------------------------------
!*Very useful safety check
!*-----------------------------------------
! PRINT*, MAXVAL(From_2Indices_To_1Index_WC)
! PRINT*, From_2Indices_To_1Index_WC
! STOP
!*-----------------------------------------



!*Neighbours only for the internal of the domain
!*->Inside the domain
indc=0
DO jj=1,nElemsY
  DO ii=1,nElemsX

    !*C,B,R,U,L
    indc=indc+1

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
    Neighbours_1Index_WC(indc,1)=indc

    !*B
    Neighbours_1Index_WC(indc,2)=From_2Indices_To_1Index_WC(ii,jj-1)

    !*R
    Neighbours_1Index_WC(indc,3)=From_2Indices_To_1Index_WC(ii+1,jj)

    !*U
    Neighbours_1Index_WC(indc,4)=From_2Indices_To_1Index_WC(ii,jj+1)

    !*L
    Neighbours_1Index_WC(indc,5)=From_2Indices_To_1Index_WC(ii-1,jj)

    !*-----------------------------------------
    !*SAFETY CHECK
    !*-----------------------------------------
    ! PRINT*
    ! PRINT*, indc 
    ! PRINT*, "C", Neighbours_1Index_WC(indc,1)
    ! PRINT*, "B", Neighbours_1Index_WC(indc,2)
    ! PRINT*, "R", Neighbours_1Index_WC(indc,3)
    ! PRINT*, "U", Neighbours_1Index_WC(indc,4)
    ! PRINT*, "L", Neighbours_1Index_WC(indc,5)
    !*-----------------------------------------

  END DO
END DO



#ifdef IMEX
#ifdef COUNTJACOBI
JacobiCounter_WC = 0
ALLOCATE(JacobiIterations_WC(1:maxLenghtJacobiCounter))
#endif

!*NB: This is strongly low order because a 5 point stencil is assumed
!*5=C,D,R,U,L
NNZsparse_WC = nElemsX*nElemsY*5 

!*We need to take BCs into account now

!------------------------------!
! Bottom Boundary Conditions   !
!------------------------------!
SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(1))
  CASE(1,2,5) ! Periodic,Transmissive,Reflecting
    NNZsparse_WC = NNZsparse_WC + 2*nElemsX
  CASE(3,4) ! Inflow,Outflow
    NNZsparse_WC = NNZsparse_WC + nElemsX
  CASE DEFAULT
    ErrorMessage = "Boundary condition B not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Right Boundary Conditions    !
!------------------------------!
SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(2))
  CASE(1,2,5) ! Periodic,Transmissive,Reflecting
    NNZsparse_WC = NNZsparse_WC + 2*nElemsY
  CASE(3,4) ! Inflow,Outflow
    NNZsparse_WC = NNZsparse_WC + nElemsY
  CASE DEFAULT
    ErrorMessage = "Boundary condition R not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Top Boundary Conditions      !
!------------------------------!
SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(3))
  CASE(1,2,5) ! Periodic,Transmissive,Reflecting
    NNZsparse_WC = NNZsparse_WC + 2*nElemsX
  CASE(3,4) ! Inflow,Outflow
    NNZsparse_WC = NNZsparse_WC + nElemsX
  CASE DEFAULT
    ErrorMessage = "Boundary condition T not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!------------------------------!
! Left Boundary Conditions     !
!------------------------------!
SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(4))
  CASE(1,2,5) ! Periodic,Transmissive,Reflecting
    NNZsparse_WC = NNZsparse_WC + 2*nElemsY
  CASE(3,4) ! Inflow,Outflow
    NNZsparse_WC = NNZsparse_WC + nElemsY
  CASE DEFAULT
    ErrorMessage = "Boundary condition L not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

ALLOCATE(Values_WC(NNZsparse_WC) )
ALLOCATE(Columns_WC(NNZsparse_WC))
ALLOCATE(RowStart_WC(NRows_WC+1) )
ALLOCATE(Diagonal_WC(NRows_WC)   )
ALLOCATE(rhs_WC(NRows_WC)        )
ALLOCATE(sol_WC(NRows_WC)        )


Values_WC   = 0.0
Columns_WC  = 0
RowStart_WC = 0
Diagonal_WC = 0.0
rhs_WC      = 0.0
sol_WC      = 0.0                                  
#endif
#endif










!*=================================================================
!* IMEX FOR ACTIVEFLUX
!*=================================================================
!*Structures for IMEX linear system
!*NB: The following is kind of analogous to what previously done
!*But it is organized slightly differently and more general

!*NB: For some nGhosts, the layers of ghost cells are nGhosts+1
!*We only care about the first layer
!*The other one will have -1 (non-physical) values in our structures


#ifdef ACTIVEFLUX
!*The following structures are essential for IMEX, 
!*but I need them also for the explicit version, for the linear system to initialize the velocity in Gresho
!*In any case, they are referred to W_X and W_Y and so they are meant to be present if ACTIVEFLUX flag is active

!*Staggering in X direction
NRows_X=(nElemsX+1)*nElemsY+2*(nElemsX+1)+2*nElemsY
ALLOCATE(From_2Indices_To_1Index_X(-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)) !*ii,jj  -> ll     !*Also for ghosts, same as W_X !*But careful with corners (-1)
ALLOCATE(From_1Index_To_2Indices_X(1:NRows_X,1:2))                                           !*ll,1:2 -> ii:jj  !*Also for ghosts excluding corners
ALLOCATE(Neighbours_1Index_X(1:((nElemsX+1)*nElemsY),1:5))                                       !*ll,1:5 ->        !*Just for internal

From_2Indices_To_1Index_X = -1 !*Initialization to -1 meaning non-physical
From_1Index_To_2Indices_X = -1 !*Initialization to -1 meaning non-physical
Neighbours_1Index_X       = -1 !*Initialization to -1 meaning non-physical

!*-------------------------
!*Numeration
!*-------------------------

indc=0
!*->Inside the domain
DO jj=1,nElemsY
  DO ii=1,nElemsX+1 
    indc=indc+1
    From_2Indices_To_1Index_X(ii,jj)=indc
    ! PRINT*, ii,jj,"->",indc
    From_1Index_To_2Indices_X(indc,1)=ii
    From_1Index_To_2Indices_X(indc,2)=jj

  END DO
END DO

!*->Down boundary
jj=0
DO ii=1,nElemsX+1 
  indc=indc+1
  From_2Indices_To_1Index_X(ii,jj)=indc
  ! PRINT*, ii,jj,"->",indc
  From_1Index_To_2Indices_X(indc,1)=ii
  From_1Index_To_2Indices_X(indc,2)=jj
END DO

!*->Right boundary
ii=nElemsX+2
DO jj=1,nElemsY
  indc=indc+1
  From_2Indices_To_1Index_X(ii,jj)=indc
  ! PRINT*, ii,jj,"->",indc
  From_1Index_To_2Indices_X(indc,1)=ii
  From_1Index_To_2Indices_X(indc,2)=jj
END DO

!*->Top boundary
jj=nElemsY+1
DO ii=1,nElemsX+1 
  indc=indc+1
  From_2Indices_To_1Index_X(ii,jj)=indc
  ! PRINT*, ii,jj,"->",indc
  From_1Index_To_2Indices_X(indc,1)=ii
  From_1Index_To_2Indices_X(indc,2)=jj
END DO

!*->Left boundary
ii=0
DO jj=1,nElemsY
  indc=indc+1
  From_2Indices_To_1Index_X(ii,jj)=indc
  ! PRINT*, ii,jj,"->",indc
  From_1Index_To_2Indices_X(indc,1)=ii
  From_1Index_To_2Indices_X(indc,2)=jj
END DO

!*-----------------------------------------
!*Very useful safety check
!*-----------------------------------------
! PRINT*, MAXVAL(From_2Indices_To_1Index_X)
! PRINT*, From_2Indices_To_1Index_X
! STOP
!*-----------------------------------------



!*Neighbours only for the internal of the domain
!*->Inside the domain
indc=0
DO jj=1,nElemsY
  DO ii=1,nElemsX+1 

    !*C,B,R,U,L
    indc=indc+1

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
    Neighbours_1Index_X(indc,1)=indc

    !*B
    Neighbours_1Index_X(indc,2)=From_2Indices_To_1Index_X(ii,jj-1)

    !*R
    Neighbours_1Index_X(indc,3)=From_2Indices_To_1Index_X(ii+1,jj)

    !*U
    Neighbours_1Index_X(indc,4)=From_2Indices_To_1Index_X(ii,jj+1)

    !*L
    Neighbours_1Index_X(indc,5)=From_2Indices_To_1Index_X(ii-1,jj)

    !*-----------------------------------------
    !*SAFETY CHECK
    !*-----------------------------------------
    ! PRINT*
    ! PRINT*, indc 
    ! PRINT*, "C", Neighbours_1Index_X(indc,1)
    ! PRINT*, "B", Neighbours_1Index_X(indc,2)
    ! PRINT*, "R", Neighbours_1Index_X(indc,3)
    ! PRINT*, "U", Neighbours_1Index_X(indc,4)
    ! PRINT*, "L", Neighbours_1Index_X(indc,5)
    !*-----------------------------------------

  END DO
END DO


!*Staggering in Y direction
NRows_Y=nElemsX*(nElemsY+1)+2*nElemsX+2*(nElemsY+1)
ALLOCATE(From_2Indices_To_1Index_Y(-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)) !*ii,jj  -> ll     !*Also for ghosts, same as W_Y !*But careful with corners (-1)
ALLOCATE(From_1Index_To_2Indices_Y(1:NRows_Y,1:2))                                           !*ll,1:2 -> ii:jj  !*Also for ghosts excluding corners
ALLOCATE(Neighbours_1Index_Y(1:(nElemsX*(nElemsY+1)),1:5))                                       !*ll,1:5 ->        !*Just for internal

From_2Indices_To_1Index_Y = -1 !*Initialization to -1 meaning non-physical
From_1Index_To_2Indices_Y = -1 !*Initialization to -1 meaning non-physical
Neighbours_1Index_Y       = -1 !*Initialization to -1 meaning non-physical

!*-------------------------
!*Numeration
!*-------------------------

indc=0
!*->Inside the domain
DO jj=1,nElemsY+1
  DO ii=1,nElemsX 
    indc=indc+1
    From_2Indices_To_1Index_Y(ii,jj)=indc
    ! PRINT*, ii,jj,"->",indc
    From_1Index_To_2Indices_Y(indc,1)=ii
    From_1Index_To_2Indices_Y(indc,2)=jj

  END DO
END DO

!*->Down boundary
jj=0
DO ii=1,nElemsX
  indc=indc+1
  From_2Indices_To_1Index_Y(ii,jj)=indc
  ! PRINT*, ii,jj,"->",indc
  From_1Index_To_2Indices_Y(indc,1)=ii
  From_1Index_To_2Indices_Y(indc,2)=jj
END DO

!*->Right boundary
ii=nElemsX+1
DO jj=1,nElemsY+1
  indc=indc+1
  From_2Indices_To_1Index_Y(ii,jj)=indc
  ! PRINT*, ii,jj,"->",indc
  From_1Index_To_2Indices_Y(indc,1)=ii
  From_1Index_To_2Indices_Y(indc,2)=jj
END DO

!*->Top boundary
jj=nElemsY+2
DO ii=1,nElemsX 
  indc=indc+1
  From_2Indices_To_1Index_Y(ii,jj)=indc
  ! PRINT*, ii,jj,"->",indc
  From_1Index_To_2Indices_Y(indc,1)=ii
  From_1Index_To_2Indices_Y(indc,2)=jj
END DO

!*->Left boundary
ii=0
DO jj=1,nElemsY+1
  indc=indc+1
  From_2Indices_To_1Index_Y(ii,jj)=indc
  ! PRINT*, ii,jj,"->",indc
  From_1Index_To_2Indices_Y(indc,1)=ii
  From_1Index_To_2Indices_Y(indc,2)=jj
END DO

!*-----------------------------------------
!*Very useful safety check
!*-----------------------------------------
! PRINT*, MAXVAL(From_2Indices_To_1Index_Y)
! PRINT*, From_2Indices_To_1Index_Y
! STOP
!*-----------------------------------------

!*Neighbours only for the internal of the domain
!*->Inside the domain
indc=0
DO jj=1,nElemsY+1
  DO ii=1,nElemsX 

    !*C,B,R,U,L
    indc=indc+1

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
    Neighbours_1Index_Y(indc,1)=indc

    !*B
    Neighbours_1Index_Y(indc,2)=From_2Indices_To_1Index_Y(ii,jj-1)

    !*R
    Neighbours_1Index_Y(indc,3)=From_2Indices_To_1Index_Y(ii+1,jj)

    !*U
    Neighbours_1Index_Y(indc,4)=From_2Indices_To_1Index_Y(ii,jj+1)

    !*L
    Neighbours_1Index_Y(indc,5)=From_2Indices_To_1Index_Y(ii-1,jj)

    !*-----------------------------------------
    !*SAFETY CHECK
    !*-----------------------------------------
    ! PRINT*
    ! PRINT*, indc 
    ! PRINT*, "C", Neighbours_1Index_Y(indc,1)
    ! PRINT*, "B", Neighbours_1Index_Y(indc,2)
    ! PRINT*, "R", Neighbours_1Index_Y(indc,3)
    ! PRINT*, "U", Neighbours_1Index_Y(indc,4)
    ! PRINT*, "L", Neighbours_1Index_Y(indc,5)
    !*-----------------------------------------

  END DO
END DO



#if defined(IMEX) || defined(IMEXMOMENTUM)
!*For staggering in X direction
#ifdef COUNTJACOBI
JacobiCounter_X = 0
ALLOCATE(JacobiIterations_X(1:maxLenghtJacobiCounter))
#endif

!*NB: This is strongly low order because a 5 point stencil is assumed
!*5=C,D,R,U,L
NNZsparse_X = (nElemsX+1)*nElemsY*5 !*NB:+1 in X direction 

!*We need to take BCs into account now

!------------------------------!
! Bottom Boundary Conditions   !
!------------------------------!
SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(1))
  CASE(1,2,5) ! Periodic,Transmissive,Reflecting
    NNZsparse_X = NNZsparse_X + 2*(nElemsX+1)
  CASE(3,4) ! Inflow,Outflow
    NNZsparse_X = NNZsparse_X + nElemsX+1
  CASE DEFAULT
    ErrorMessage = "Boundary condition B not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Right Boundary Conditions    !
!------------------------------!
SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(2))
  CASE(1,2,5) ! Periodic,Transmissive,Reflecting
    NNZsparse_X = NNZsparse_X + 2*nElemsY
  CASE(3,4) ! Inflow,Outflow
    NNZsparse_X = NNZsparse_X + nElemsY
  CASE DEFAULT
    ErrorMessage = "Boundary condition R not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Top Boundary Conditions      !
!------------------------------!
SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(3))
  CASE(1,2,5) ! Periodic,Transmissive,Reflecting
    NNZsparse_X = NNZsparse_X + 2*(nElemsX+1)
  CASE(3,4) ! Inflow,Outflow
    NNZsparse_X = NNZsparse_X + nElemsX+1
  CASE DEFAULT
    ErrorMessage = "Boundary condition T not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!------------------------------!
! Left Boundary Conditions     !
!------------------------------!
SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(4))
  CASE(1,2,5) ! Periodic,Transmissive,Reflecting
    NNZsparse_X = NNZsparse_X + 2*nElemsY
  CASE(3,4) ! Inflow,Outflow
    NNZsparse_X = NNZsparse_X + nElemsY
  CASE DEFAULT
    ErrorMessage = "Boundary condition L not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

ALLOCATE(Values_X(NNZsparse_X) )
ALLOCATE(Columns_X(NNZsparse_X))
ALLOCATE(RowStart_X(NRows_X+1) )
ALLOCATE(Diagonal_X(NRows_X)   )
ALLOCATE(rhs_X(NRows_X)        )
ALLOCATE(sol_X(NRows_X)        )


Values_X   = 0.0
Columns_X  = 0
RowStart_X = 0
Diagonal_X = 0.0
rhs_X      = 0.0
sol_X      = 0.0

!*For staggering in Y direction
#ifdef COUNTJACOBI
JacobiCounter_Y = 0
ALLOCATE(JacobiIterations_Y(1:maxLenghtJacobiCounter))
#endif

!*NB: This is strongly low order because a 5 point stencil is assumed
!*5=C,D,R,U,L
NNZsparse_Y = nElemsX*(nElemsY+1)*5 !*NB:+1 in Y direction 

!*We need to take BCs into account now

!------------------------------!
! Bottom Boundary Conditions   !
!------------------------------!
SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(1))
  CASE(1,2,5) ! Periodic,Transmissive,Reflecting
    NNZsparse_Y = NNZsparse_Y + 2*nElemsX
  CASE(3,4) ! Inflow,Outflow
    NNZsparse_Y = NNZsparse_Y + nElemsX
  CASE DEFAULT
    ErrorMessage = "Boundary condition B not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Right Boundary Conditions    !
!------------------------------!
SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(2))
  CASE(1,2,5) ! Periodic,Transmissive,Reflecting
    NNZsparse_Y = NNZsparse_Y + 2*(nElemsY+1)
  CASE(3,4) ! Inflow,Outflow
    NNZsparse_Y = NNZsparse_Y + nElemsY+1
  CASE DEFAULT
    ErrorMessage = "Boundary condition R not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Top Boundary Conditions      !
!------------------------------!
SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(3))
  CASE(1,2,5) ! Periodic,Transmissive,Reflecting
    NNZsparse_Y = NNZsparse_Y + 2*nElemsX
  CASE(3,4) ! Inflow,Outflow
    NNZsparse_Y = NNZsparse_Y + nElemsX
  CASE DEFAULT
    ErrorMessage = "Boundary condition T not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!------------------------------!
! Left Boundary Conditions     !
!------------------------------!
SELECT CASE(BoundaryConditionsTypeIMEXLinearSystem(4))
  CASE(1,2,5) ! Periodic,Transmissive,Reflecting
    NNZsparse_Y = NNZsparse_Y + 2*(nElemsY+1)
  CASE(3,4) ! Inflow,Outflow
    NNZsparse_Y = NNZsparse_Y + nElemsY+1
  CASE DEFAULT
    ErrorMessage = "Boundary condition L not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


ALLOCATE(Values_Y(NNZsparse_Y) )
ALLOCATE(Columns_Y(NNZsparse_Y))
ALLOCATE(RowStart_Y(NRows_Y+1) )
ALLOCATE(Diagonal_Y(NRows_Y)   )
ALLOCATE(rhs_Y(NRows_Y)        )
ALLOCATE(sol_Y(NRows_Y)        )

Values_Y   = 0.0
Columns_Y  = 0
RowStart_Y = 0
Diagonal_Y = 0.0
rhs_Y      = 0.0
sol_Y      = 0.0
                                  
#endif
#endif



!-------------------------------------------------------------------------------!
END SUBROUTINE InitializeFiniteVolume
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE FillInitialConditions()
!-------------------------------------------------------------------------------!
USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Equation,           ONLY: ExactFunction
USE MOD_Equation,           ONLY: Gravitational_Potential   
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Gravitational_Potential_Averages
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGP
USE MOD_FiniteVolume2D_vars,ONLY: MESH_SX
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
#endif

#ifdef ACTIVEFLUX
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_X

USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_Y

USE MOD_Equation           ,ONLY: BoundaryConditions
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(nVar, nGPs, nGPs) :: Utemp, Wtemp
REAL, DIMENSION(nGPs, nGPs)       :: GravPottemp
INTEGER :: ii, jj, iGP, jGP, iVar
REAL, DIMENSION(nVar)             :: debugu,debugv
REAL, DIMENSION(0:nElemsX+1)      :: xcoord
REAL, DIMENSION(0:nElemsY+1)      :: ycoord
REAL                              ::rho1, uu1, vv1, pp1
REAL                              ::rho2, uu2, vv2, pp2
REAL                              :: dx,dxx,dy,dyy
INTEGER                           :: i,j,k
!-------------------------------------------------------------------------------!

U     = 0.
V     = 0.
Utemp = 0.
GravPottemp = 0.

#ifdef CENTEREDPRIMITIVE
WC = 0.0
#endif

DO jj=1,nElemsY
  DO ii=1,nElemsX
    DO iGP=1,nGPs
      DO jGP=1,nGPs
        ! compute cell average of conservative variables
        CALL ExactFunction(&
          InitialCondition,0.0,MeshGP(:,ii,jj,iGP,jGP),Utemp(1:nVar,iGP,jGP))

        U(1:nVar, ii, jj) = U(1:nVar, ii, jj) + WeightsGP(iGP,jGP) * Utemp(1:nVar,iGP,jGP)

#ifdef CENTEREDPRIMITIVE
        CALL ConsToPrim(Utemp(1:nVar,iGP,jGP),Wtemp(1:nVar,iGP,jGP))
        WC(1:nVar, ii, jj) = WC(1:nVar, ii, jj) + WeightsGP(iGP,jGP) * Wtemp(1:nVar,iGP,jGP)
#endif


        ! compute cell average of bathymetry            
        GravPottemp(iGP,jGP) = Gravitational_Potential(MeshGP(:,ii,jj,iGP,jGP))
        Gravitational_Potential_Averages(ii,jj)    = Gravitational_Potential_Averages(ii,jj) + WeightsGP(iGP,jGP) * GravPottemp(iGP,jGP)
      END DO
    END DO

    CALL ConsToPrim(U(1:nVar,ii,jj),V(1:nVar,ii,jj))
  END DO
END DO

#if(1==1)
!*====================
!*DIRECT INTIALIZATION FOR IMPLOSION TEST TO OVERCOME MP EFFECTS
!*====================
IF (InitialCondition .EQ. 10) THEN
      dx  =MESH_SX(1)/REAL(nElemsX)
      dy  =MESH_SX(2)/REAL(nElemsY)
      dxx =0.50*dx
      dyy =0.50*dy

      do j=0,nElemsX+1
       xcoord(j)=MESH_X0(1)+REAL(j)*dx-dxx
      enddo
      do k=0,nElemsY+1
       ycoord(k)=MESH_X0(2)+REAL(k)*dy-dyy
      enddo
      rho1 = 0.1250
      uu1  = 0.0
      vv1  = 0.0
      pp1  = 0.140

      rho2 = 1.0
      uu2  = 0.0
      vv2  = 0.0
      pp2  = 1.0

      do k=1,nElemsY
       do j=1,nElemsX
        if((abs(xcoord(j))+abs(ycoord(k)))<0.150)then
         u(1,j,k)=rho1
         u(2,j,k)=rho1*uu1
         u(3,j,k)=rho1*vv1
         u(4,j,k)=pp1/(Gmm-1.0)+0.50*rho1*(uu1*uu1+vv1*vv1)
        else
         u(1,j,k)=rho2
         u(2,j,k)=rho2*uu2
         u(3,j,k)=rho2*vv2
         u(4,j,k)=pp2/(Gmm-1.0)+0.50*rho2*(uu2*uu2+vv2*vv2)
        endif
#ifdef CENTEREDPRIMITIVE
        CALL ConsToPrim(U(1:nVar,j,k),WC(1:nVar,j,k))
#endif
       enddo
      enddo

END IF
!*====================
#endif


#ifdef ACTIVEFLUX
!*Compute cell averages of primitive variables on staggered meshes

!*For staggering in X direction
W_X=0.
Wtemp=0.

DO jj=1,nElemsY
  DO ii=1,nElemsX+1 !*NB: +1 in X direction
    DO iGP=1,nGPs
      DO jGP=1,nGPs
        CALL ExactFunction(&
          InitialCondition,0.0,MeshGP_X(:,ii,jj,iGP,jGP),Utemp(1:nVar,iGP,jGP))

        CALL ConsToPrim(Utemp(1:nVar,iGP,jGP),Wtemp(1:nVar,iGP,jGP))

        W_X(1:nVar, ii, jj) = W_X(1:nVar, ii, jj) + WeightsGP(iGP,jGP) * Wtemp(1:nVar,iGP,jGP)

        !*-------------------------------------
        !*SAFETY CHECK
        !*-------------------------------------
        ! CALL ExactFunction(&
        !   InitialCondition,0.0,MeshGP_X(:,ii,jj,iGP,jGP),debugu)
        ! CALL ConsToPrim(debugu,debugv)
        ! DO iVar=1,nVar
        !   IF( debugv(iVar) .NE. Wtemp(iVar,iGP,jGP) ) THEN
        !     PRINT*, debugv(iVar) - Wtemp(iVar,iGP,jGP)
        !     PRINT*, "Initialization problem in Mesh staggered X"
        !     STOP
        !   END IF
        ! END DO
        ! PRINT*, jj, ii, W_X(:,ii,jj) !*,MeshGP_X(:,ii,jj,iGP,jGP)
        !*-------------------------------------

      END DO
    END DO

  END DO
END DO

!*For staggering in Y direction
W_Y=0.
Wtemp=0.

DO jj=1,nElemsY+1 !*NB: +1 in Y direction
  DO ii=1,nElemsX
    DO iGP=1,nGPs
      DO jGP=1,nGPs
        CALL ExactFunction(&
          InitialCondition,0.0,MeshGP_Y(:,ii,jj,iGP,jGP),Utemp(1:nVar,iGP,jGP))

        CALL ConsToPrim(Utemp(1:nVar,iGP,jGP),Wtemp(1:nVar,iGP,jGP))

        W_Y(1:nVar, ii, jj) = W_Y(1:nVar, ii, jj) + WeightsGP(iGP,jGP) * Wtemp(1:nVar,iGP,jGP)

        !*-------------------------------------
        !*SAFETY CHECK
        !*-------------------------------------
        ! CALL ExactFunction(&
        !   InitialCondition,0.0,MeshGP_Y(:,ii,jj,iGP,jGP),debugu)
        ! CALL ConsToPrim(debugu,debugv)
        ! DO iVar=1,nVar
        !   IF( debugv(iVar) .NE. Wtemp(iVar,iGP,jGP) ) THEN
        !     PRINT*, debugv(iVar) - Wtemp(iVar,iGP,jGP)
        !     PRINT*, "Initialization problem in Mesh staggered Y"
        !     STOP
        !   END IF
        ! END DO
        ! PRINT*, jj, ii, W_Y(:,ii,jj) !*, MeshGP_Y(:,ii,jj,iGP,jGP)
        !*-------------------------------------

      END DO
    END DO

  END DO
END DO

#if(1==0)
#ifdef RELAXATION
PRINT*, "In the end divergence free preparation of data was not helping. I'm stopping here."
STOP
IF (InitialCondition .EQ. 100) THEN
  !*For Gresho vortex, divergence free preparation of data
  CALL BoundaryConditions(0.0)
  CALL Divergence_Free_Preparation_Of_Data_W_X()
  CALL Divergence_Free_Preparation_Of_Data_W_Y()
  CALL BoundaryConditions(0.0)
END IF
#endif
#endif

#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE FillInitialConditions
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE InitializeWBVariables()
!-------------------------------------------------------------------------------!
USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: UtWB
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: FXWB
USE MOD_FiniteVolume2D_vars,ONLY: FX
USE MOD_FiniteVolume2D_vars,ONLY: FYWB
USE MOD_FiniteVolume2D_vars,ONLY: FY
USE MOD_FiniteVolume2D_vars,ONLY: SWB
USE MOD_FiniteVolume2D_vars,ONLY: S
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGP
USE MOD_Equation           ,ONLY: ExactFunctionWB
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(nVar, nGPs, nGPs) :: Utemp
INTEGER :: ii, jj, iGP, jGP
!-------------------------------------------------------------------------------!

! compute UWB
Utemp = 0.
DO jj=1,nElemsY
  DO ii=1,nElemsX
    DO iGP=1,nGPs
      DO jGP=1,nGPs
        CALL ExactFunctionWB(InitialCondition,MeshGP(:,ii,jj,iGP,jGP),Utemp(1:nVar,iGP,jGP))
        U(1:nVar, ii, jj) = U(1:nVar, ii, jj) + WeightsGP(iGP,jGP) * Utemp(1:nVar,iGP,jGP)
      END DO
    END DO
    CALL ConsToPrim(U(1:nVar,ii,jj),V(1:nVar,ii,jj))
  END DO
END DO

CALL FVTimeDerivative(0.)

UtWB = Ut
FXWB = FX
FYWB = FY
SWB  = S


!-------------------------------------------------------------------------------!
END SUBROUTINE InitializeWBVariables
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE FVTimeDerivative(t,Optional_Input)
!-------------------------------------------------------------------------------!
USE MOD_Equation,       ONLY: BoundaryConditions
USE MOD_Equation,       ONLY: SourceTerms
USE MOD_Equation,       ONLY: SourceTerm_Conserved_INPUT_CONSERVED
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
USE MOD_Equation,       ONLY: SourceTerm_Primitive_INPUT_PRIMITIVE
#endif
USE MOD_Reconstruction, ONLY: ReconstructionX
USE MOD_Reconstruction, ONLY: ReconstructionY
USE MOD_Reconstruction, ONLY: ReconstructionFixX
USE MOD_Reconstruction, ONLY: ReconstructionFixY
USE MOD_ShocksIndicator,ONLY: ShocksIndicatorX
USE MOD_ShocksIndicator,ONLY: ShocksIndicatorY
#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars, ONLY: WC
USE MOD_FiniteVolume2D_vars, ONLY: WCt
#ifdef PATHCONSERVATIVESHOCKDETECTION
USE MOD_FiniteVolume2D_vars, ONLY: Cell_To_Limit_PCSD
#endif
#endif
#ifdef ACTIVEFLUX
USE MOD_Reconstruction, ONLY: FromCellAveragesToPointValues

!*Unified subroutines
USE MOD_FiniteVolume2D_vars, ONLY: W_X
USE MOD_FiniteVolume2D_vars, ONLY: W_Y
USE MOD_FiniteVolume2D_vars, ONLY: Wt_X
USE MOD_FiniteVolume2D_vars, ONLY: Wt_Y
#endif
USE MOD_Reconstruction,      ONLY: ReconstructionX_INPUT_PRIMITIVE
USE MOD_Reconstruction,      ONLY: ReconstructionY_INPUT_PRIMITIVE
USE MOD_Reconstruction,      ONLY: ReconstructionX_INPUT_CONSERVED
USE MOD_Reconstruction,      ONLY: ReconstructionY_INPUT_CONSERVED
USE MOD_Reconstruction,      ONLY: ReconstructionFixX_INPUT_PRIMITIVE
USE MOD_Reconstruction,      ONLY: ReconstructionFixY_INPUT_PRIMITIVE
USE MOD_Reconstruction,      ONLY: ReconstructionFixX_INPUT_CONSERVED
USE MOD_Reconstruction,      ONLY: ReconstructionFixY_INPUT_CONSERVED
USE MOD_FiniteVolume2D_vars, ONLY: nGhosts
USE MOD_FiniteVolume2D_vars, ONLY: nElemsX
USE MOD_FiniteVolume2D_vars, ONLY: nElemsY
USE MOD_FiniteVolume2D_vars, ONLY: U
USE MOD_FiniteVolume2D_vars, ONLY: nVar
#ifdef MULTIFLUID
USE MOD_FiniteVolume2D_vars, ONLY: Troubled_Cell_W_X
USE MOD_FiniteVolume2D_vars, ONLY: Troubled_Cell_W_Y
#endif
USE MOD_FiniteVolume2D_vars, ONLY: nGhosts
USE MOD_FiniteVolume2D_vars, ONLY: nGPs
USE MOD_FiniteVolume2D_vars, ONLY: WM
USE MOD_FiniteVolume2D_vars, ONLY: WP
USE MOD_FiniteVolume2D_vars, ONLY: MeshBary
USE MOD_FiniteVolume2D_vars, ONLY: MeshGP
USE MOD_FiniteVolume2D_vars, ONLY: nDims
#ifdef ACTIVEFLUX
USE MOD_FiniteVolume2D_vars, ONLY: MeshBary_X
USE MOD_FiniteVolume2D_vars, ONLY: MeshGP_X
USE MOD_FiniteVolume2D_vars, ONLY: MeshBary_Y
USE MOD_FiniteVolume2D_vars, ONLY: MeshGP_Y
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,    INTENT(IN)           :: t
INTEGER, INTENT(IN), OPTIONAL :: Optional_Input
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
#ifdef CENTEREDPRIMITIVE
REAL            :: W_support( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Exactly as U
#endif
#ifdef ACTIVEFLUX
REAL            :: W_X_support( 1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)!*NB:+1 in X direction with respect to U
REAL            :: W_Y_support( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)!*NB:+1 in Y direction with respect to U
#endif
!*For passing from conserved to primitive
REAL            :: WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1) 
REAL            :: WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1) 
INTEGER         :: ii,jj


CALL BoundaryConditions(t)




#if defined(CENTEREDPRIMITIVE)
!*=========================
!*U
!*=========================
#ifdef DETECTEACHSTAGE
IF (.NOT.(PRESENT(Optional_Input))) THEN
  CALL Shock_Detector_Based_On_Path_Conservative()
  ! PRINT*, "Detect"
ELSE
  ! PRINT*, "Not detect"
END IF
#endif


#ifdef RECONSTRUCTFROMCONSERVED
CALL ComputeTransitionMatricesConservedInput(U( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),  &  !*Vector
                          & -nGhosts,nElemsX+nGhosts+1, &  !*X bounds vector  
                          & -nGhosts,nElemsY+nGhosts+1,   &  !*Y bounds vector
                          & 1,nElemsX, &                   !*X bounds loop    
                          & 1,nElemsY)                       !*Y bounds loop 
#else
CALL ComputeTransitionMatricesPrimitiveInput(WC( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),  &  !*Vector
                          & -nGhosts,nElemsX+nGhosts+1, &  !*X bounds vector  
                          & -nGhosts,nElemsY+nGhosts+1,   &  !*Y bounds vector
                          & 1,nElemsX, &                   !*X bounds loop    
                          & 1,nElemsY)                       !*Y bounds loop 
#endif

#ifdef RECONSTRUCTFROMCONSERVED
CALL ReconstructionX_INPUT_CONSERVED(U( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&    !*Y bounds vector  
                          & 0,nElemsX+1, &                 !*X bounds loop  
                          & 1,nElemsY)                   !*Y bounds loop    

#ifdef PATHCONSERVATIVESHOCKDETECTION
CALL ReconstructionFixX_INPUT_CONSERVED(U( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&    !*Y bounds vector  
                          & 0,nElemsX+1, &                 !*X bounds loop  
                          & 1,nElemsY, &                   !*Y bounds loop    
                          & Cell_To_Limit_PCSD( -nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) )
#endif

CALL ConsToPrimInTheWholeMeshInterfaces(WM(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1),0,nElemsX+1+1,0,nElemsY+1+1,WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1))
WM(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)=WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)
CALL ConsToPrimInTheWholeMeshInterfaces(WP(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1),0,nElemsX+1+1,0,nElemsY+1+1,WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1))
WP(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)=WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)

#else
CALL ReconstructionX_INPUT_PRIMITIVE(WC( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&    !*Y bounds vector  
                          & 0,nElemsX+1, &                 !*X bounds loop  
                          & 1,nElemsY)                   !*Y bounds loop    

#ifdef PATHCONSERVATIVESHOCKDETECTION
CALL ReconstructionFixX_INPUT_PRIMITIVE(WC( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&    !*Y bounds vector  
                          & 0,nElemsX+1, &                 !*X bounds loop  
                          & 1,nElemsY, &                   !*Y bounds loop    
                          & Cell_To_Limit_PCSD( -nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) )
#endif

#endif

CALL NumericalFluxFX()

#ifdef RECONSTRUCTFROMCONSERVED
CALL ReconstructionY_INPUT_CONSERVED(U( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&  !*Y bounds vector 
                          & 1,nElemsX, &                   !*X bounds loop
                          & 0,nElemsY+1)                 !*Y bounds loop    


#ifdef PATHCONSERVATIVESHOCKDETECTION
CALL ReconstructionFixY_INPUT_CONSERVED(U( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&  !*Y bounds vector 
                          & 1,nElemsX, &                   !*X bounds loop
                          & 0,nElemsY+1, &                 !*Y bounds loop    
                          & Cell_To_Limit_PCSD( -nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) )
#endif


CALL ConsToPrimInTheWholeMeshInterfaces(WM(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1),0,nElemsX+1+1,0,nElemsY+1+1,WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1))
WM(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)=WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)
CALL ConsToPrimInTheWholeMeshInterfaces(WP(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1),0,nElemsX+1+1,0,nElemsY+1+1,WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1))
WP(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)=WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)


#else
CALL ReconstructionY_INPUT_PRIMITIVE(WC( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&  !*Y bounds vector 
                          & 1,nElemsX, &                   !*X bounds loop
                          & 0,nElemsY+1)                 !*Y bounds loop    

#ifdef PATHCONSERVATIVESHOCKDETECTION
CALL ReconstructionFixY_INPUT_PRIMITIVE(WC( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&  !*Y bounds vector 
                          & 1,nElemsX, &                   !*X bounds loop
                          & 0,nElemsY+1, &                 !*Y bounds loop    
                          & Cell_To_Limit_PCSD( -nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) )
#endif

#endif



CALL NumericalFluxFY()

CALL SourceTerm_Conserved_INPUT_CONSERVED(t,U(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&    !*Y bounds vector  
                          & 1,nElemsX, &                   !*X bounds loop  
                          & 1,nElemsY, &                   !*Y bounds loop    
                          & MeshBary(1:nDims,1:nElemsX,1:nElemsY), &
                          & 1,nElemsX, &
                          & 1,nElemsY, &
                          & MeshGP   (1:nDims,1:nElemsX,1:nElemsY,1:nGPs,1:nGPs), &
                          & 1,nElemsX, &
                          & 1,nElemsY)

CALL UpdateTimeDerivative()

!*=========================
!*WC
!*=========================
#ifndef FIRSTORDERPRIMITIVE

#ifdef RECONSTRUCTFROMCONSERVED
CALL ReconstructionX_INPUT_CONSERVED(U( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&    !*Y bounds vector  
                          & 0,nElemsX+1, &                 !*X bounds loop  
                          & 1,nElemsY)                   !*Y bounds loop    

#ifdef PATHCONSERVATIVESHOCKDETECTION
CALL ReconstructionFixX_INPUT_CONSERVED(U( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&    !*Y bounds vector  
                          & 0,nElemsX+1, &                 !*X bounds loop  
                          & 1,nElemsY, &                   !*Y bounds loop    
                          & Cell_To_Limit_PCSD( -nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) )
#endif

CALL ConsToPrimInTheWholeMeshInterfaces(WM(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1),0,nElemsX+1+1,0,nElemsY+1+1,WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1))
WM(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)=WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)
CALL ConsToPrimInTheWholeMeshInterfaces(WP(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1),0,nElemsX+1+1,0,nElemsY+1+1,WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1))
WP(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)=WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)
#else
CALL ReconstructionX_INPUT_PRIMITIVE(WC( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&    !*Y bounds vector  
                          & 0,nElemsX+1, &                 !*X bounds loop  
                          & 1,nElemsY)                   !*Y bounds loop    

#ifdef PATHCONSERVATIVESHOCKDETECTION
CALL ReconstructionFixX_INPUT_PRIMITIVE(WC( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&    !*Y bounds vector  
                          & 0,nElemsX+1, &                 !*X bounds loop  
                          & 1,nElemsY, &                   !*Y bounds loop    
                          & Cell_To_Limit_PCSD( -nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) )
#endif

#endif


CALL NumericalFluxFX_W_INPUT( 0,nElemsX,&                    !*X bounds loop    
                            & 1,nElemsY)                     !*Y bounds loop 

CALL BSurfX_W_INPUT( 0,nElemsX,&                             !*X bounds loop 
                   & 1,nElemsY)                              !*Y bounds loop 

CALL BVolX_W_INPUT(WC( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                                    !*Vector
                  & -nGhosts,nElemsX+nGhosts+1, &            !*X bounds vector  
                  & -nGhosts,nElemsY+nGhosts+1,   &          !*Y bounds vector
                  & 1,nElemsX, &                             !*X bounds loop    
                  & 1,nElemsY)                               !*Y bounds loop  

#ifdef RECONSTRUCTFROMCONSERVED
CALL ReconstructionY_INPUT_CONSERVED(U( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&  !*Y bounds vector 
                          & 1,nElemsX, &                   !*X bounds loop
                          & 0,nElemsY+1)                 !*Y bounds loop    

#ifdef PATHCONSERVATIVESHOCKDETECTION
CALL ReconstructionFixY_INPUT_CONSERVED(U( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&  !*Y bounds vector 
                          & 1,nElemsX, &                   !*X bounds loop
                          & 0,nElemsY+1, &                 !*Y bounds loop    
                          & Cell_To_Limit_PCSD( -nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) )
#endif


CALL ConsToPrimInTheWholeMeshInterfaces(WM(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1),0,nElemsX+1+1,0,nElemsY+1+1,WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1))
WM(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)=WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)
CALL ConsToPrimInTheWholeMeshInterfaces(WP(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1),0,nElemsX+1+1,0,nElemsY+1+1,WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1))
WP(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)=WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)
#else
CALL ReconstructionY_INPUT_PRIMITIVE(WC( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&  !*Y bounds vector 
                          & 1,nElemsX, &                   !*X bounds loop
                          & 0,nElemsY+1)                 !*Y bounds loop    

#ifdef PATHCONSERVATIVESHOCKDETECTION
CALL ReconstructionFixY_INPUT_PRIMITIVE(WC( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&  !*Y bounds vector 
                          & 1,nElemsX, &                   !*X bounds loop
                          & 0,nElemsY+1, &                 !*Y bounds loop    
                          & Cell_To_Limit_PCSD( -nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) )
#endif

#endif

CALL NumericalFluxFY_W_INPUT( 1,nElemsX, &                   !*X bounds loop    
                            & 0,nElemsY )                    !*Y bounds loop  

CALL BSurfY_W_INPUT( 1,nElemsX,&                             !*X bounds loop 
                   & 0,nElemsY)                              !*Y bounds loop 

CALL BVolY_W_INPUT(WC( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                                    !*Vector
                  & -nGhosts,nElemsX+nGhosts+1, &            !*X bounds vector  
                  & -nGhosts,nElemsY+nGhosts+1,   &          !*Y bounds vector
                  & 1,nElemsX, &                             !*X bounds loop    
                  & 1,nElemsY)                               !*Y bounds loop  


CALL SourceTerm_Primitive_INPUT_PRIMITIVE(t,WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&    !*Y bounds vector  
                          & 1,nElemsX, &                   !*X bounds loop  
                          & 1,nElemsY, &                   !*Y bounds loop    
                          & MeshBary(1:nDims,1:nElemsX,1:nElemsY), &
                          & 1,nElemsX, &
                          & 1,nElemsY, &
                          & MeshGP   (1:nDims,1:nElemsX,1:nElemsY,1:nGPs,1:nGPs), &
                          & 1,nElemsX, &
                          & 1,nElemsY)


CALL UpdateTimeDerivative_W_INPUT(WCt(1:nVar,1:nElemsX,1:nElemsY), &    !*Vector
                                 & 1,nElemsX, &                         !*X bounds loop    
                                 & 1,nElemsY)                           !*Y bounds loop  

#endif

#elif defined(ACTIVEFLUX)

#ifdef MULTIFLUID
CALL Identify_Interface_Cells()
#endif


W_X_support=W_X
W_Y_support=W_Y
#ifdef MOMENTUMINPRIMITIVEVARIABLES
#ifdef RECONSTRUCTIONVELOCITY
DO jj=-nGhosts,nElemsY+nGhosts+1
  DO ii=-nGhosts,nElemsX+nGhosts+1+1
    IF (W_X_support(1,ii,jj) .NE. 0.0) THEN
      W_X_support(2:3,ii,jj)=W_X_support(2:3,ii,jj)/W_X_support(1,ii,jj)
    END IF
  END DO
END DO
DO jj=-nGhosts,nElemsY+nGhosts+1+1
  DO ii=-nGhosts,nElemsX+nGhosts+1
    IF (W_Y_support(1,ii,jj) .NE. 0.0) THEN
      W_Y_support(2:3,ii,jj)=W_Y_support(2:3,ii,jj)/W_Y_support(1,ii,jj)
    END IF
  END DO
END DO
#endif
#endif

#ifdef MULTIFLUID
#ifdef RECONSTRUCTIONPHI
DO jj=-nGhosts,nElemsY+nGhosts+1
  DO ii=-nGhosts,nElemsX+nGhosts+1+1
    IF (W_X_support(1,ii,jj) .NE. 0.0) THEN
      W_X_support(nVar,ii,jj)=W_X_support(nVar,ii,jj)/W_X_support(1,ii,jj)
    END IF
  END DO
END DO
DO jj=-nGhosts,nElemsY+nGhosts+1+1
  DO ii=-nGhosts,nElemsX+nGhosts+1
    IF (W_Y_support(1,ii,jj) .NE. 0.0) THEN
      W_Y_support(nVar,ii,jj)=W_Y_support(nVar,ii,jj)/W_Y_support(1,ii,jj)
    END IF
  END DO
END DO
#endif
#endif


!*-------------------------
!*U
!*-------------------------
!*From cell averages W_X, W_Y to W_qp
!*FluxFX
!*FluxFY
!*SourceTerms
!*UpdateTimeDerivative(U)
!*-------------------------
#ifndef PRIMITIVEONLY

#if(1==1)
CALL FromCellAveragesToPointValues()
CALL FluxFX()                !*Filling FX for the update of U through reconstructed values
CALL FluxFY()                !*Filling FY for the update of U through reconstructed values
CALL SourceTerms(t)          !*Filling S for the update of U. Standard computation.
CALL UpdateTimeDerivative()  !*Filling Ut
#else
PRINT*, "Just for debugging, proper evolution of conserved variables"
STOP
CALL ComputeTransitionMatricesConservedInput(U( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),  &  !*Vector
                          & -nGhosts,nElemsX+nGhosts+1, &  !*X bounds vector  
                          & -nGhosts,nElemsY+nGhosts+1,   &  !*Y bounds vector
                          & 1,nElemsX, &                   !*X bounds loop    
                          & 1,nElemsY)                       !*Y bounds loop 

CALL ShocksIndicatorX()
CALL ReconstructionX()
CALL ReconstructionFixX()
CALL PositivityLimiterX()

CALL ConsToPrimInTheWholeMeshInterfaces(WM(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1),0,nElemsX+1+1,0,nElemsY+1+1,WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1))
WM(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)=WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)
CALL ConsToPrimInTheWholeMeshInterfaces(WP(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1),0,nElemsX+1+1,0,nElemsY+1+1,WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1))
WP(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)=WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)


CALL NumericalFluxFX()

CALL ShocksIndicatorY() 
CALL ReconstructionY()
CALL ReconstructionFixY()
CALL PositivityLimiterY()

CALL ConsToPrimInTheWholeMeshInterfaces(WM(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1),0,nElemsX+1+1,0,nElemsY+1+1,WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1))
WM(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)=WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)
CALL ConsToPrimInTheWholeMeshInterfaces(WP(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1),0,nElemsX+1+1,0,nElemsY+1+1,WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1))
WP(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)=WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)

CALL NumericalFluxFY()

CALL SourceTerms(t)
CALL UpdateTimeDerivative()
#endif

#endif

!*-------------------------
!*W_X
!*RMK: W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
!*-------------------------
!*ReconstructionX_W_X
!*NumericalFluxFX_W_X
!*BSurfX_W_X
!*BVolX_W_X
!*---
!*ReconstructionY_W_X
!*NumericalFluxFY_W_X
!*BSurfY_W_X
!*BVolY_W_X
!*---
!*UpdateTimeDerivative_W_X
!*-------------------------
CALL ComputeTransitionMatricesPrimitiveInput(W_X( 1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),  &  !*Vector
                          & -nGhosts,nElemsX+nGhosts+1+1, &  !*X bounds vector  !*NB:+1 in X direction
                          & -nGhosts,nElemsY+nGhosts+1,   &  !*Y bounds vector
                          & 1,nElemsX+1, &                   !*X bounds loop    !*NB:+1 in X direction
                          & 1,nElemsY)                       !*Y bounds loop 


#ifdef PWCINTROUBLEDCELLS
#ifdef MULTIFLUID
CALL ReconstructionX_INPUT_PRIMITIVE(W_X_support( 1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1), &                            !*Vector
                          & -nGhosts,nElemsX+nGhosts+1+1, &  !*X bounds vector  !*NB:+1 in X direction
                          & -nGhosts,nElemsY+nGhosts+1,   &  !*Y bounds vector
                          & 0,nElemsX+1+1, &                 !*X bounds loop    !*NB:+1 in X direction
                          & 1,nElemsY, &                     !*Y bounds loop  
                          & Troubled_Cell_W_X)
#else
PRINT*, "Impossible to use PWCINTROUBLEDCELLS"
PRINT*, "Disable this flag"
PRINT*, "You should not be here"
STOP
#endif
#else
CALL ReconstructionX_INPUT_PRIMITIVE(W_X_support( 1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1), &                            !*Vector
                          & -nGhosts,nElemsX+nGhosts+1+1, &  !*X bounds vector  !*NB:+1 in X direction
                          & -nGhosts,nElemsY+nGhosts+1,   &  !*Y bounds vector
                          & 0,nElemsX+1+1, &                 !*X bounds loop    !*NB:+1 in X direction
                          & 1,nElemsY)                       !*Y bounds loop  
#endif


CALL Adjust_Interface_values()


CALL NumericalFluxFX_W_INPUT( 0,nElemsX+1,&                  !*X bounds loop    !*NB:+1 in X direction
                            & 1,nElemsY)                     !*Y bounds loop 


CALL BSurfX_W_INPUT( 0,nElemsX+1,&                           !*X bounds loop !*NB:+1 in X direction, like NumericalFluxFX_W_X
                   & 1,nElemsY)                              !*Y bounds loop 

#ifdef PWCINTROUBLEDCELLS
#ifdef MULTIFLUID
CALL BVolX_W_INPUT(W_X( 1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1), &                                    !*Vector
                  & -nGhosts,nElemsX+nGhosts+1+1, &          !*X bounds vector  !*NB:+1 in X direction
                  & -nGhosts,nElemsY+nGhosts+1,   &          !*Y bounds vector
                  & 1,nElemsX+1, &                           !*X bounds loop    !*NB:+1 in X direction
                  & 1,nElemsY, &                             !*Y bounds loop  
                  & Troubled_Cell_W_X)
#else
PRINT*, "Impossible to use PWCINTROUBLEDCELLS"
PRINT*, "Disable this flag"
PRINT*, "You should not be here"
STOP
#endif
#else
CALL BVolX_W_INPUT(W_X( 1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1), &                                    !*Vector
                  & -nGhosts,nElemsX+nGhosts+1+1, &          !*X bounds vector  !*NB:+1 in X direction
                  & -nGhosts,nElemsY+nGhosts+1,   &          !*Y bounds vector
                  & 1,nElemsX+1, &                           !*X bounds loop    !*NB:+1 in X direction
                  & 1,nElemsY)                               !*Y bounds loop  
#endif

#ifdef PWCINTROUBLEDCELLS
#ifdef MULTIFLUID
CALL ReconstructionY_INPUT_PRIMITIVE(W_X_support( 1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1), &                            !*Vector
                          & -nGhosts,nElemsX+nGhosts+1+1, &  !*X bounds vector  !*NB:+1 in X direction
                          & -nGhosts,nElemsY+nGhosts+1,   &  !*Y bounds vector
                          & 1,nElemsX+1, &                   !*X bounds loop    !*NB:+1 in X direction
                          & 0,nElemsY+1, &                   !*Y bounds loop  
                          & Troubled_Cell_W_X)
#else
PRINT*, "Impossible to use PWCINTROUBLEDCELLS"
PRINT*, "Disable this flag"
PRINT*, "You should not be here"
STOP
#endif
#else
CALL ReconstructionY_INPUT_PRIMITIVE(W_X_support( 1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1), &                            !*Vector
                          & -nGhosts,nElemsX+nGhosts+1+1, &  !*X bounds vector  !*NB:+1 in X direction
                          & -nGhosts,nElemsY+nGhosts+1,   &  !*Y bounds vector
                          & 1,nElemsX+1, &                   !*X bounds loop    !*NB:+1 in X direction
                          & 0,nElemsY+1)                     !*Y bounds loop  
#endif


CALL Adjust_Interface_values()

CALL NumericalFluxFY_W_INPUT( 1,nElemsX+1, &                 !*X bounds loop    !*NB:+1 in X direction
                            & 0,nElemsY )                    !*Y bounds loop  

CALL BSurfY_W_INPUT( 1,nElemsX+1,&                           !*X bounds loop !*NB:+1 in X direction
                   & 0,nElemsY)                              !*Y bounds loop 

#ifdef PWCINTROUBLEDCELLS
#ifdef MULTIFLUID
CALL BVolY_W_INPUT(W_X( 1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1), &                                    !*Vector
                  & -nGhosts,nElemsX+nGhosts+1+1, &          !*X bounds vector  !*NB:+1 in X direction
                  & -nGhosts,nElemsY+nGhosts+1,   &          !*Y bounds vector
                  & 1,nElemsX+1, &                           !*X bounds loop    !*NB:+1 in X direction
                  & 1,nElemsY, &                             !*Y bounds loop  
                  & Troubled_Cell_W_X)
#else
PRINT*, "Impossible to use PWCINTROUBLEDCELLS"
PRINT*, "Disable this flag"
PRINT*, "You should not be here"
STOP
#endif
#else
CALL BVolY_W_INPUT(W_X( 1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1), &                                    !*Vector
                  & -nGhosts,nElemsX+nGhosts+1+1, &          !*X bounds vector  !*NB:+1 in X direction
                  & -nGhosts,nElemsY+nGhosts+1,   &          !*Y bounds vector
                  & 1,nElemsX+1, &                           !*X bounds loop    !*NB:+1 in X direction
                  & 1,nElemsY)                               !*Y bounds loop  
#endif

CALL SourceTerm_Primitive_INPUT_PRIMITIVE(t,W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1+1,  &  !*X bounds vector !*NB:+1 in X direction
                          &-nGhosts,nElemsY+nGhosts+1,&    !*Y bounds vector  
                          & 1,nElemsX+1, &                   !*X bounds loop !*NB:+1 in X direction
                          & 1,nElemsY, &                   !*Y bounds loop    
                          & MeshBary_X(1:nDims,1:nElemsX+1,1:nElemsY), &
                          & 1,nElemsX+1, &                  !*NB:+1 in X direction
                          & 1,nElemsY, &
                          & MeshGP_X(1:nDims,1:nElemsX+1,1:nElemsY,1:nGPs,1:nGPs), &
                          & 1,nElemsX+1, &                   !*NB:+1 in X direction
                          & 1,nElemsY)

!*Filling Wt_X
CALL UpdateTimeDerivative_W_INPUT(Wt_X(1:nVar,1:nElemsX+1,1:nElemsY), &   !*Vector
                                 & 1,nElemsX+1, &                         !*X bounds loop    !*NB:+1 in X direction
                                 & 1,nElemsY)                             !*Y bounds loop  
                                  

!*-------------------------
!*W_Y
!*RMK: W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
!*-------------------------
!*ReconstructionX_W_Y
!*NumericalFluxFX_W_Y
!*BSurfX_W_Y
!*BVolX_W_Y
!*---
!*ReconstructionY_W_Y
!*NumericalFluxFY_W_Y
!*BSurfY_W_Y
!*BVolY_W_Y
!*---
!*UpdateTimeDerivative_W_Y
!*-------------------------

CALL ComputeTransitionMatricesPrimitiveInput(W_Y( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1+1,&  !*Y bounds vector  !*NB:+1 in Y direction
                          & 1,nElemsX, &                 !*X bounds loop  
                          & 1,nElemsY+1)                 !*Y bounds loop    !*NB:+1 in Y direction


#ifdef PWCINTROUBLEDCELLS
#ifdef MULTIFLUID
CALL ReconstructionX_INPUT_PRIMITIVE(W_Y_support( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1+1,&  !*Y bounds vector  !*NB:+1 in Y direction
                          & 0,nElemsX+1, &                 !*X bounds loop  
                          & 1,nElemsY+1, &                 !*Y bounds loop    !*NB:+1 in Y direction
                          & Troubled_Cell_W_Y)
#else
PRINT*, "Impossible to use PWCINTROUBLEDCELLS"
PRINT*, "Disable this flag"
PRINT*, "You should not be here"
STOP
#endif
#else
CALL ReconstructionX_INPUT_PRIMITIVE(W_Y_support( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1+1,&  !*Y bounds vector  !*NB:+1 in Y direction
                          & 0,nElemsX+1, &                 !*X bounds loop  
                          & 1,nElemsY+1)                   !*Y bounds loop    !*NB:+1 in Y direction
#endif



CALL Adjust_Interface_values()


CALL NumericalFluxFX_W_INPUT( 0,nElemsX,&                  !*X bounds loop    
                            & 1,nElemsY+1)                 !*Y bounds loop    !*NB:+1 in Y direction


CALL BSurfX_W_INPUT( 0,nElemsX,&                           !*X bounds loop 
                   & 1,nElemsY+1)                          !*Y bounds loop !*NB:+1 in Y direction, like NumericalFluxFX_W_Y


#ifdef PWCINTROUBLEDCELLS
#ifdef MULTIFLUID
CALL BVolX_W_INPUT(W_Y( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1), &                                  !*Vector
                  &-nGhosts,nElemsX+nGhosts+1,  &          !*X bounds vector
                  &-nGhosts,nElemsY+nGhosts+1+1,&          !*Y bounds vector  !*NB:+1 in Y direction
                  & 1,nElemsX, &                           !*X bounds loop  
                  & 1,nElemsY+1, &                         !*Y bounds loop    !*NB:+1 in Y direction
                  & Troubled_Cell_W_Y)
#else
PRINT*, "Impossible to use PWCINTROUBLEDCELLS"
PRINT*, "Disable this flag"
PRINT*, "You should not be here"
STOP
#endif
#else
CALL BVolX_W_INPUT(W_Y( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1), &                                  !*Vector
                  &-nGhosts,nElemsX+nGhosts+1,  &          !*X bounds vector
                  &-nGhosts,nElemsY+nGhosts+1+1,&          !*Y bounds vector  !*NB:+1 in Y direction
                  & 1,nElemsX, &                           !*X bounds loop  
                  & 1,nElemsY+1)                           !*Y bounds loop    !*NB:+1 in Y direction
#endif


#ifdef PWCINTROUBLEDCELLS
#ifdef MULTIFLUID
CALL ReconstructionY_INPUT_PRIMITIVE(W_Y_support( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1+1,&  !*Y bounds vector  !*NB:+1 in Y direction
                          & 1,nElemsX, &                   !*X bounds loop
                          & 0,nElemsY+1+1, &               !*Y bounds loop    !*NB:+1 in Y direction
                          & Troubled_Cell_W_Y)
#else
PRINT*, "Impossible to use PWCINTROUBLEDCELLS"
PRINT*, "Disable this flag"
PRINT*, "You should not be here"
STOP
#endif
#else
CALL ReconstructionY_INPUT_PRIMITIVE(W_Y_support( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1+1,&  !*Y bounds vector  !*NB:+1 in Y direction
                          & 1,nElemsX, &                   !*X bounds loop
                          & 0,nElemsY+1+1)                 !*Y bounds loop    !*NB:+1 in Y direction
#endif


CALL Adjust_Interface_values()


CALL NumericalFluxFY_W_INPUT( 1,nElemsX, &                 !*X bounds loop
                            & 0,nElemsY+1 )                !*Y bounds loop    !*NB:+1 in Y direction  


CALL BSurfY_W_INPUT( 1,nElemsX, &                          !*X bounds loop
                   & 0,nElemsY+1)                          !*Y bounds loop    !*NB:+1 in Y direction

#ifdef PWCINTROUBLEDCELLS
#ifdef MULTIFLUID
CALL BVolY_W_INPUT(W_Y( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1), &                                  !*Vector
                  &-nGhosts,nElemsX+nGhosts+1,  &          !*X bounds vector
                  &-nGhosts,nElemsY+nGhosts+1+1,&          !*Y bounds vector  !*NB:+1 in Y direction
                  & 1,nElemsX, &                           !*X bounds loop  
                  & 1,nElemsY+1, &                         !*Y bounds loop    !*NB:+1 in Y direction
                  & Troubled_Cell_W_Y)
#else
PRINT*, "Impossible to use PWCINTROUBLEDCELLS"
PRINT*, "Disable this flag"
PRINT*, "You should not be here"
STOP
#endif
#else
CALL BVolY_W_INPUT(W_Y( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1), &                                  !*Vector
                  &-nGhosts,nElemsX+nGhosts+1,  &          !*X bounds vector
                  &-nGhosts,nElemsY+nGhosts+1+1,&          !*Y bounds vector  !*NB:+1 in Y direction
                  & 1,nElemsX, &                           !*X bounds loop  
                  & 1,nElemsY+1)                           !*Y bounds loop    !*NB:+1 in Y direction
#endif

CALL SourceTerm_Primitive_INPUT_PRIMITIVE(t,W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1+1,&    !*Y bounds vector !*NB:+1 in Y direction  
                          & 1,nElemsX, &                   !*X bounds loop
                          & 1,nElemsY+1, &                   !*Y bounds loop !*NB:+1 in Y direction    
                          & MeshBary_Y(1:nDims,1:nElemsX,1:nElemsY+1), &
                          & 1,nElemsX, &                  
                          & 1,nElemsY+1, &                  !*NB:+1 in Y direction
                          & MeshGP_Y(1:nDims,1:nElemsX,1:nElemsY+1,1:nGPs,1:nGPs), &
                          & 1,nElemsX, &                   
                          & 1,nElemsY+1)                      !*NB:+1 in Y direction


!*<========  !*Filling Wt_Y
CALL UpdateTimeDerivative_W_INPUT(Wt_Y(1:nVar,1:nElemsX,1:nElemsY+1), &  !*Vector
                                 & 1,nElemsX, &                          !*X bounds loop
                                 & 1,nElemsY+1)                          !*Y bounds loop    !*NB:+1 in Y direction  


#else
!*-------------------------
!*RMK: U(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
!*-------------------------
CALL ComputeTransitionMatricesConservedInput(U( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),  &  !*Vector
                          & -nGhosts,nElemsX+nGhosts+1, &  !*X bounds vector  
                          & -nGhosts,nElemsY+nGhosts+1,   &  !*Y bounds vector
                          & 1,nElemsX, &                   !*X bounds loop    
                          & 1,nElemsY)                       !*Y bounds loop 

#if(1==1)
CALL ReconstructionX_INPUT_CONSERVED(U( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&    !*Y bounds vector  
                          & 0,nElemsX+1, &                 !*X bounds loop  
                          & 1,nElemsY)                   !*Y bounds loop    
#else
CALL ShocksIndicatorX()
CALL ReconstructionX()
CALL ReconstructionFixX()
CALL PositivityLimiterX()
#endif
CALL ConsToPrimInTheWholeMeshInterfaces(WM(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1),0,nElemsX+1+1,0,nElemsY+1+1,WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1))
WM(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)=WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)
CALL ConsToPrimInTheWholeMeshInterfaces(WP(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1),0,nElemsX+1+1,0,nElemsY+1+1,WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1))
WP(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)=WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)


CALL NumericalFluxFX()

#if(1==1)
CALL ReconstructionY_INPUT_CONSERVED(U( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&  !*Y bounds vector 
                          & 1,nElemsX, &                   !*X bounds loop
                          & 0,nElemsY+1)                 !*Y bounds loop    
#else
CALL ShocksIndicatorY() 
CALL ReconstructionY()
CALL ReconstructionFixY()
CALL PositivityLimiterY()
#endif
CALL ConsToPrimInTheWholeMeshInterfaces(WM(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1),0,nElemsX+1+1,0,nElemsY+1+1,WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1))
WM(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)=WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)
CALL ConsToPrimInTheWholeMeshInterfaces(WP(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1),0,nElemsX+1+1,0,nElemsY+1+1,WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1))
WP(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)=WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)


CALL NumericalFluxFY()

CALL SourceTerms(t)
CALL UpdateTimeDerivative()

#endif


!-------------------------------------------------------------------------------!
END SUBROUTINE FVTimeDerivative
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE UpdateTimeDerivative()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: Ut, S
USE MOD_Mesh,               ONLY: GlobalElem
USE MOD_FiniteVolume2D_vars,ONLY: UtWB
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: FX
USE MOD_FiniteVolume2D_vars,ONLY: FY
USE MOD_FiniteVolume2D_vars,ONLY: FXWB
USE MOD_FiniteVolume2D_vars,ONLY: FYWB
USE MOD_FiniteVolume2D_vars,ONLY: SWB

#ifdef PATANKAR
USE MOD_FiniteVolume2D_vars,ONLY: ProductionSparse
USE MOD_FiniteVolume2D_vars,ONLY: DestructionSparse
USE MOD_FiniteVolume2D_vars,ONLY: ColumnsVector
USE MOD_FiniteVolume2D_vars,ONLY: RowStart
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, K, L, iSparseVect, iRowStart
!-------------------------------------------------------------------------------!

!Production = 0.
!Destruction = 0.

#ifdef PATANKAR
  iSparseVect = 0
  iRowStart=0
#endif  

#ifdef WELLBALANCED
  FX=FX-FXWB
  FY=FY-FYWB
  S =S -SWB
#endif

DO jj=1,nElemsY
  DO ii=1,nElemsX
#ifdef PATANKAR   
    K = GlobalElem(ii,jj)
    
    iSparseVect = iSparseVect + 1
    iRowStart = iRowStart +1

    L = GlobalElem(ii,jj-1)
    ProductionSparse(iSparseVect) = max(FY(1,ii+0,jj-1),0.)/Mesh_DX(2)
    DestructionSparse(iSparseVect) = -min(FY(1,ii+0,jj-1),0.)/Mesh_DX(2)

    iSparseVect = iSparseVect + 1
    L = GlobalElem(ii-1,jj)
    ProductionSparse(iSparseVect) =   max(FX(1,ii-1,jj+0),0.)/Mesh_DX(1)
    DestructionSparse(iSparseVect) = -min(FX(1,ii-1,jj+0),0.)/Mesh_DX(1)


    iSparseVect = iSparseVect + 1
    L = GlobalElem(ii,jj)
    ProductionSparse(iSparseVect) = 0.
    DestructionSparse(iSparseVect) = 0.


    iSparseVect = iSparseVect + 1
    L = GlobalElem(ii+1,jj)
    ProductionSparse(iSparseVect) =  max(-FX(1,ii+0,jj+0),0.)/Mesh_DX(1)
    DestructionSparse(iSparseVect) = -min(-FX(1,ii+0,jj+0),0.)/Mesh_DX(1)

    iSparseVect = iSparseVect + 1
    L = GlobalElem(ii,jj+1)
    ProductionSparse(iSparseVect) = max(-FY(1,ii+0,jj+0),0.)/Mesh_DX(2)
    DestructionSparse(iSparseVect) = -min(-FY(1,ii+0,jj+0),0.)/Mesh_DX(2)

    Ut(2:nVar,ii,jj) = S(2:nVar,ii,jj) &
                     - (FX(2:nVar,ii+0,jj+0)-FX(2:nVar,ii-1,jj+0))/Mesh_DX(1) &
                     - (FY(2:nVar,ii+0,jj+0)-FY(2:nVar,ii+0,jj-1))/Mesh_DX(2)

#else
    Ut(1:nVar,ii,jj) = S(1:nVar,ii,jj) &
                     - (FX(1:nVar,ii+0,jj+0)-FX(1:nVar,ii-1,jj+0))/Mesh_DX(1) &
                     - (FY(1:nVar,ii+0,jj+0)-FY(1:nVar,ii+0,jj-1))/Mesh_DX(2)

#endif
  END DO !ii
END DO !jj


!-------------------------------------------------------------------------------!
END SUBROUTINE UpdateTimeDerivative
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE NumericalFluxFX()
!-------------------------------------------------------------------------------!
USE MOD_Equation,           ONLY: RiemannSolver
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: FX
USE MOD_FiniteVolume2D_vars,ONLY: FluxX
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: NormVectX, TangVectX
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGPBnd
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, iGP
!-------------------------------------------------------------------------------!

FX    = 0.0
FluxX = 0.0

DO jj=1,nElemsY
  DO ii=0,nElemsX
    CALL RiemannSolver(&
                       WP(1:nVar,1:nGPs,ii+0,jj),&
                       WM(1:nVar,1:nGPs,ii+1,jj),&
                       NormVectX(1:nDims),&
                       TangVectX(1:nDims),&
                       FluxX(1:nVar,1:nGPs,ii,jj))
  END DO
END DO

DO iGP = 1,nGPs
  FX(:,:,:) = FX(:,:,:) + weightsGPBnd(iGP)*FluxX(:,iGP,:,:)
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE NumericalFluxFX
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE NumericalFluxFY()
!-------------------------------------------------------------------------------!
USE MOD_Equation,           ONLY: RiemannSolver
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: FY
USE MOD_FiniteVolume2D_vars,ONLY: FluxY
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: NormVectY, TangVectY
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGPBnd
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, iGP
!-------------------------------------------------------------------------------!

FY    = 0.0
FluxY = 0.0

DO jj=0,nElemsY
  DO ii=1,nElemsX
    CALL RiemannSolver(&
                       WP(1:nVar,1:nGPs,ii,jj+0),&
                       WM(1:nVar,1:nGPs,ii,jj+1),&
                       NormVectY(1:nDims),&
                       TangVectY(1:nDims),&
                       FluxY(1:nVar,1:nGPs,ii,jj))
  END DO
END DO

DO iGP = 1,nGPs
  FY(:,:,:) = FY(:,:,:) + weightsGPBnd(iGP)*FluxY(:,iGP,:,:)
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE NumericalFluxFY
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE PositivityLimiterX()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGPBnd
USE MOD_FiniteVolume2D_vars,ONLY: wLobatto    
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DENSITY   
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, jGP
REAL    :: alpha, beta, csiX, theta, mmin 
!-------------------------------------------------------------------------------!


DO jj=1,nElemsY
  DO ii=0,nElemsX+1
     alpha = 0.  
     beta  = 0.
     DO jGP = 1,nGPs
        alpha = alpha + WeightsGPBnd(jGP)*WM(1,jGP,ii+0,jj) 
        beta  = beta  + WeightsGPBnd(jGP)*WP(1,jGP,ii+0,jj)
     END DO
     csiX = ( U(1,ii,jj) - wLobatto*alpha - wLobatto*beta )/( 1. - 2.*wLobatto ) 
     DO jGP = 1,nGPs
        mmin = MIN( csiX , WM(1,jGP,ii+0,jj), WP(1,jGP,ii+0,jj))
        IF ( U(1,ii,jj) .EQ. mmin ) THEN
          theta = 1. 
        ELSE
          theta = MIN( 1. , ABS( (U(1,ii,jj)-MIN_DENSITY)/(U(1,ii,jj)-mmin) ) )
        ENDIF
        WP(1,jGP,ii+0,jj) = U(1,ii,jj) + theta * ( WP(1,jGP,ii+0,jj) - U(1,ii,jj) ) 
        WM(1,jGP,ii+0,jj) = U(1,ii,jj) + theta * ( WM(1,jGP,ii+0,jj) - U(1,ii,jj) )
     END DO
  END DO
END DO



END SUBROUTINE PositivityLimiterX
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE PositivityLimiterY()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: FY
USE MOD_FiniteVolume2D_vars,ONLY: FluxY
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: NormVectY, TangVectY
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGPBnd
USE MOD_FiniteVolume2D_vars,ONLY: wLobatto    
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DENSITY   
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, jGP
REAL    :: alpha, beta, csiY, theta, mmin 
!-------------------------------------------------------------------------------!


DO jj=0,nElemsY+1
  DO ii=1,nElemsX
     alpha = 0.  
     beta  = 0.
     DO jGP = 1,nGPs
        alpha = alpha + WeightsGPBnd(jGP)*WM(1,jGP,ii,jj+0) 
        beta  = beta  + WeightsGPBnd(jGP)*WP(1,jGP,ii,jj+0)
     END DO
     csiY = ( U(1,ii,jj) - wLobatto*alpha - wLobatto*beta )/( 1. - 2.*wLobatto ) 
     DO jGP = 1,nGPs
        mmin = MIN( csiY , WM(1,jGP,ii,jj+0), WP(1,jGP,ii,jj+0))
        IF ( U(1,ii,jj) .EQ. mmin ) THEN
          theta = 1. 
        ELSE
          theta = MIN( 1. , ABS( (U(1,ii,jj)-MIN_DENSITY)/(U(1,ii,jj)-mmin) ) )
        ENDIF
        WP(1,jGP,ii,jj+0) = U(1,ii,jj) + theta * ( WP(1,jGP,ii,jj+0) - U(1,ii,jj) ) 
        WM(1,jGP,ii,jj+0) = U(1,ii,jj) + theta * ( WM(1,jGP,ii,jj+0) - U(1,ii,jj) )
     END DO
  END DO
END DO


END SUBROUTINE PositivityLimiterY
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE PrimToConsInTheWholeMeshCells(W,vbXs,vbXe,vbYs,vbYe,Wout)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_Equation,           ONLY: PrimToCons
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,vbXs:vbXe,vbYs:vbYe),    INTENT(IN)          :: W
INTEGER,                                        INTENT(IN)          :: vbXs,vbXe,vbYs,vbYe !*Bounds for the vector
REAL, DIMENSION(1:nVar,vbXs:vbXe,vbYs:vbYe),    INTENT(INOUt)       :: Wout
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, jGP
REAL    :: alpha, beta, csiY, theta, mmin 
!-------------------------------------------------------------------------------!

DO jj=vbYs,vbYe
  DO ii=vbXs,vbXe
    IF (W(1,ii,jj) .NE. 0.0) THEN !*If density is not equal to 0.0
      CALL PrimToCons(W(1:nVar,ii,jj),Wout(1:nVar,ii,jj))
    END IF
  END DO
END DO


END SUBROUTINE PrimToConsInTheWholeMeshCells
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ConsToPrimInTheWholeMeshCells(W,vbXs,vbXe,vbYs,vbYe,Wout)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_Equation,           ONLY: ConsToPrim
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,vbXs:vbXe,vbYs:vbYe),    INTENT(IN)          :: W
INTEGER,                                        INTENT(IN)          :: vbXs,vbXe,vbYs,vbYe !*Bounds for the vector
REAL, DIMENSION(1:nVar,vbXs:vbXe,vbYs:vbYe),    INTENT(INOUt)       :: Wout
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, jGP
REAL    :: alpha, beta, csiY, theta, mmin 
!-------------------------------------------------------------------------------!

DO jj=vbYs,vbYe
  DO ii=vbXs,vbXe
    IF (W(1,ii,jj) .NE. 0.0) THEN !*If density is not equal to 0.0
      CALL ConsToPrim(W(1:nVar,ii,jj),Wout(1:nVar,ii,jj))
    END IF
  END DO
END DO


END SUBROUTINE ConsToPrimInTheWholeMeshCells
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE PrimToConsInTheWholeMeshInterfaces(W,vbXs,vbXe,vbYs,vbYe,Wout)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_Equation,           ONLY: PrimToCons
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,1:nGPs,vbXs:vbXe,vbYs:vbYe),    INTENT(IN)          :: W
INTEGER,                                               INTENT(IN)          :: vbXs,vbXe,vbYs,vbYe !*Bounds for the vector
REAL, DIMENSION(1:nVar,1:nGPs,vbXs:vbXe,vbYs:vbYe),    INTENT(INOUt)       :: Wout
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, jGP
REAL    :: alpha, beta, csiY, theta, mmin 
!-------------------------------------------------------------------------------!

DO jj=vbYs,vbYe
  DO ii=vbXs,vbXe
    DO jGP=1,nGPs
      IF (W(1,jGP,ii,jj) .NE. 0.0) THEN !*If density is not equal to 0.0
        CALL PrimToCons(W(1:nVar,jGP,ii,jj),Wout(1:nVar,jGP,ii,jj))
      END IF
    END DO
  END DO
END DO


END SUBROUTINE PrimToConsInTheWholeMeshInterfaces
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ConsToPrimInTheWholeMeshInterfaces(W,vbXs,vbXe,vbYs,vbYe,Wout)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_Equation,           ONLY: ConsToPrim
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,1:nGPs,vbXs:vbXe,vbYs:vbYe),    INTENT(IN)          :: W
INTEGER,                                               INTENT(IN)          :: vbXs,vbXe,vbYs,vbYe !*Bounds for the vector
REAL, DIMENSION(1:nVar,1:nGPs,vbXs:vbXe,vbYs:vbYe),    INTENT(INOUT)       :: Wout
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, jGP
REAL    :: alpha, beta, csiY, theta, mmin 
!-------------------------------------------------------------------------------!

DO jj=vbYs,vbYe
  DO ii=vbXs,vbXe
    DO jGP=1,nGPs
      IF (W(1,jGP,ii,jj) .NE. 0.0) THEN !*If density is not equal to 0.0
        CALL ConsToPrim(W(1:nVar,jGP,ii,jj),Wout(1:nVar,jGP,ii,jj))
      END IF
    END DO
  END DO
END DO


END SUBROUTINE ConsToPrimInTheWholeMeshInterfaces
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ComputeTransitionMatricesPrimitiveInput(W,vbXs,vbXe,vbYs,vbYe,lbXs,lbXe,lbYs,lbYe)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: LLX
USE MOD_FiniteVolume2D_vars,ONLY: RRX
USE MOD_FiniteVolume2D_vars,ONLY: LLY
USE MOD_FiniteVolume2D_vars,ONLY: RRY
USE MOD_FiniteVolume2D_vars,ONLY: MCtoP
USE MOD_FiniteVolume2D_vars,ONLY: MPtoC
USE MOD_FiniteVolume2D_vars,ONLY: NormVectX
USE MOD_FiniteVolume2D_vars,ONLY: NormVectY
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX, nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_Equation,           ONLY: LeftEigenvectorsXDirectionPrimitiveSystem
USE MOD_Equation,           ONLY: RightEigenvectorsXDirectionPrimitiveSystem
USE MOD_Equation,           ONLY: LeftEigenvectorsYDirectionPrimitiveSystem
USE MOD_Equation,           ONLY: RightEigenvectorsYDirectionPrimitiveSystem
USE MOD_Equation,           ONLY: LeftEigenvectorsRotationalInvariancePrimitiveSystem
USE MOD_Equation,           ONLY: RightEigenvectorsRotationalInvariancePrimitiveSystem
USE MOD_Equation,           ONLY: TransitionMatrixConsToPrim
USE MOD_Equation,           ONLY: TransitionMatrixPrimToCons
USE MOD_Equation,           ONLY: PrimToCons
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,vbXs:vbXe,vbYs:vbYe),    INTENT(IN)          :: W                   !*NB: Conserved variables assumed as inputs
INTEGER,                                        INTENT(IN)          :: vbXs,vbXe,vbYs,vbYe !*Bounds for the vector
INTEGER,                                        INTENT(IN)          :: lbXs,lbXe,lbYs,lbYe !*Bounds for the loop   
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, indi,indj
REAL    :: LMat(1:nVar,1:nVar)
REAL    :: RMat(1:nVar,1:nVar)
REAL    :: MMat(1:nVar,1:nVar)
REAL    :: Cons(1:nVar)
!-------------------------------------------------------------------------------!

#if(1==1)

  DO ii=lbXs,lbXe
    DO jj=lbYs,lbYe
      ! PRINT*, ii, jj, W(1:nVar,ii,jj)
      CALL PrimToCons(W(1:nVar,ii,jj),Cons(1:nVar))
      CALL LeftEigenvectorsXDirectionPrimitiveSystem(          W(1:nVar,ii,jj),LLX(:,:,ii,jj))
      CALL RightEigenvectorsXDirectionPrimitiveSystem(         W(1:nVar,ii,jj),RRX(:,:,ii,jj))
      CALL LeftEigenvectorsYDirectionPrimitiveSystem(          W(1:nVar,ii,jj),LLY(:,:,ii,jj))
      CALL RightEigenvectorsYDirectionPrimitiveSystem(         W(1:nVar,ii,jj),RRY(:,:,ii,jj))
      CALL TransitionMatrixConsToPrim(Cons(1:nVar),MCtoP(:,:,ii,jj))
      CALL TransitionMatrixPrimToCons(Cons(1:nVar),MPtoC(:,:,ii,jj))
    END DO
  END DO 

  !*Ghost cells (only one layer)
  !*B,R,U,L

  !*B
  jj=lbYs-1
  DO ii=lbXs,lbXe
      ! PRINT*, ii, jj, "B", W(1:nVar,ii,jj)
      CALL PrimToCons(W(1:nVar,ii,jj),Cons(1:nVar))
      CALL LeftEigenvectorsXDirectionPrimitiveSystem(          W(1:nVar,ii,jj),LLX(:,:,ii,jj))
      CALL RightEigenvectorsXDirectionPrimitiveSystem(         W(1:nVar,ii,jj),RRX(:,:,ii,jj))
      CALL LeftEigenvectorsYDirectionPrimitiveSystem(          W(1:nVar,ii,jj),LLY(:,:,ii,jj))
      CALL RightEigenvectorsYDirectionPrimitiveSystem(         W(1:nVar,ii,jj),RRY(:,:,ii,jj))
      CALL TransitionMatrixConsToPrim(Cons(1:nVar),MCtoP(:,:,ii,jj))
      CALL TransitionMatrixPrimToCons(Cons(1:nVar),MPtoC(:,:,ii,jj))
  END DO 

  !*R
  ii=lbXe+1
  DO jj=lbYs,lbYe
      ! PRINT*, ii, jj, "R", W(1:nVar,ii,jj)
      CALL PrimToCons(W(1:nVar,ii,jj),Cons(1:nVar))
      CALL LeftEigenvectorsXDirectionPrimitiveSystem(          W(1:nVar,ii,jj),LLX(:,:,ii,jj))
      CALL RightEigenvectorsXDirectionPrimitiveSystem(         W(1:nVar,ii,jj),RRX(:,:,ii,jj))
      CALL LeftEigenvectorsYDirectionPrimitiveSystem(          W(1:nVar,ii,jj),LLY(:,:,ii,jj))
      CALL RightEigenvectorsYDirectionPrimitiveSystem(         W(1:nVar,ii,jj),RRY(:,:,ii,jj))
      CALL TransitionMatrixConsToPrim(Cons(1:nVar),MCtoP(:,:,ii,jj))
      CALL TransitionMatrixPrimToCons(Cons(1:nVar),MPtoC(:,:,ii,jj))
  END DO

  !*U
  jj=lbYe+1
  DO ii=lbXs,lbXe
      ! PRINT*, ii, jj, "U", W(1:nVar,ii,jj)
      CALL PrimToCons(W(1:nVar,ii,jj),Cons(1:nVar))
      CALL LeftEigenvectorsXDirectionPrimitiveSystem(          W(1:nVar,ii,jj),LLX(:,:,ii,jj))
      CALL RightEigenvectorsXDirectionPrimitiveSystem(         W(1:nVar,ii,jj),RRX(:,:,ii,jj))
      CALL LeftEigenvectorsYDirectionPrimitiveSystem(          W(1:nVar,ii,jj),LLY(:,:,ii,jj))
      CALL RightEigenvectorsYDirectionPrimitiveSystem(         W(1:nVar,ii,jj),RRY(:,:,ii,jj))
      CALL TransitionMatrixConsToPrim(Cons(1:nVar),MCtoP(:,:,ii,jj))
      CALL TransitionMatrixPrimToCons(Cons(1:nVar),MPtoC(:,:,ii,jj))
  END DO 

  !*L
  ii=lbXs-1
  DO jj=lbYs,lbYe
      ! PRINT*, ii, jj, "L", W(1:nVar,ii,jj)
      CALL PrimToCons(W(1:nVar,ii,jj),Cons(1:nVar))
      CALL LeftEigenvectorsXDirectionPrimitiveSystem(          W(1:nVar,ii,jj),LLX(:,:,ii,jj))
      CALL RightEigenvectorsXDirectionPrimitiveSystem(         W(1:nVar,ii,jj),RRX(:,:,ii,jj))
      CALL LeftEigenvectorsYDirectionPrimitiveSystem(          W(1:nVar,ii,jj),LLY(:,:,ii,jj))
      CALL RightEigenvectorsYDirectionPrimitiveSystem(         W(1:nVar,ii,jj),RRY(:,:,ii,jj))
      CALL TransitionMatrixConsToPrim(Cons(1:nVar),MCtoP(:,:,ii,jj))
      CALL TransitionMatrixPrimToCons(Cons(1:nVar),MPtoC(:,:,ii,jj))
  END DO

#else


  DO ii=lbXs,lbXe
    DO jj=lbYs,lbYe
      ! PRINT*, ii, jj, W(1:nVar,ii,jj)
      CALL PrimToCons(W(1:nVar,ii,jj),Cons(1:nVar))
      CALL LeftEigenvectorsRotationalInvariancePrimitiveSystem(          W(1:nVar,ii,jj),NormVectX,LLX(:,:,ii,jj))
      CALL RightEigenvectorsRotationalInvariancePrimitiveSystem(         W(1:nVar,ii,jj),NormVectX,RRX(:,:,ii,jj))
      CALL LeftEigenvectorsRotationalInvariancePrimitiveSystem(          W(1:nVar,ii,jj),NormVectY,LLY(:,:,ii,jj))
      CALL RightEigenvectorsRotationalInvariancePrimitiveSystem(         W(1:nVar,ii,jj),NormVectY,RRY(:,:,ii,jj))
      CALL TransitionMatrixConsToPrim(Cons(1:nVar),MCtoP(:,:,ii,jj))
      CALL TransitionMatrixPrimToCons(Cons(1:nVar),MPtoC(:,:,ii,jj))
    END DO
  END DO 

  !*Ghost cells (only one layer)
  !*B,R,U,L

  !*B
  jj=lbYs-1
  DO ii=lbXs,lbXe
      ! PRINT*, ii, jj, "B", W(1:nVar,ii,jj)
      CALL PrimToCons(W(1:nVar,ii,jj),Cons(1:nVar))
      CALL LeftEigenvectorsRotationalInvariancePrimitiveSystem(          W(1:nVar,ii,jj),NormVectX,LLX(:,:,ii,jj))
      CALL RightEigenvectorsRotationalInvariancePrimitiveSystem(         W(1:nVar,ii,jj),NormVectX,RRX(:,:,ii,jj))
      CALL LeftEigenvectorsRotationalInvariancePrimitiveSystem(          W(1:nVar,ii,jj),NormVectY,LLY(:,:,ii,jj))
      CALL RightEigenvectorsRotationalInvariancePrimitiveSystem(         W(1:nVar,ii,jj),NormVectY,RRY(:,:,ii,jj))
      CALL TransitionMatrixConsToPrim(Cons(1:nVar),MCtoP(:,:,ii,jj))
      CALL TransitionMatrixPrimToCons(Cons(1:nVar),MPtoC(:,:,ii,jj))
  END DO 

  !*R
  ii=lbXe+1
  DO jj=lbYs,lbYe
      ! PRINT*, ii, jj, "R", W(1:nVar,ii,jj)
      CALL PrimToCons(W(1:nVar,ii,jj),Cons(1:nVar))
      CALL LeftEigenvectorsRotationalInvariancePrimitiveSystem(          W(1:nVar,ii,jj),NormVectX,LLX(:,:,ii,jj))
      CALL RightEigenvectorsRotationalInvariancePrimitiveSystem(         W(1:nVar,ii,jj),NormVectX,RRX(:,:,ii,jj))
      CALL LeftEigenvectorsRotationalInvariancePrimitiveSystem(          W(1:nVar,ii,jj),NormVectY,LLY(:,:,ii,jj))
      CALL RightEigenvectorsRotationalInvariancePrimitiveSystem(         W(1:nVar,ii,jj),NormVectY,RRY(:,:,ii,jj))
      CALL TransitionMatrixConsToPrim(Cons(1:nVar),MCtoP(:,:,ii,jj))
      CALL TransitionMatrixPrimToCons(Cons(1:nVar),MPtoC(:,:,ii,jj))
  END DO

  !*U
  jj=lbYe+1
  DO ii=lbXs,lbXe
      ! PRINT*, ii, jj, "U", W(1:nVar,ii,jj)
      CALL PrimToCons(W(1:nVar,ii,jj),Cons(1:nVar))
      CALL LeftEigenvectorsRotationalInvariancePrimitiveSystem(          W(1:nVar,ii,jj),NormVectX,LLX(:,:,ii,jj))
      CALL RightEigenvectorsRotationalInvariancePrimitiveSystem(         W(1:nVar,ii,jj),NormVectX,RRX(:,:,ii,jj))
      CALL LeftEigenvectorsRotationalInvariancePrimitiveSystem(          W(1:nVar,ii,jj),NormVectY,LLY(:,:,ii,jj))
      CALL RightEigenvectorsRotationalInvariancePrimitiveSystem(         W(1:nVar,ii,jj),NormVectY,RRY(:,:,ii,jj))
      CALL TransitionMatrixConsToPrim(Cons(1:nVar),MCtoP(:,:,ii,jj))
      CALL TransitionMatrixPrimToCons(Cons(1:nVar),MPtoC(:,:,ii,jj))
  END DO 

  !*L
  ii=lbXs-1
  DO jj=lbYs,lbYe
      ! PRINT*, ii, jj, "L", W(1:nVar,ii,jj)
      CALL PrimToCons(W(1:nVar,ii,jj),Cons(1:nVar))
      CALL LeftEigenvectorsRotationalInvariancePrimitiveSystem(          W(1:nVar,ii,jj),NormVectX,LLX(:,:,ii,jj))
      CALL RightEigenvectorsRotationalInvariancePrimitiveSystem(         W(1:nVar,ii,jj),NormVectX,RRX(:,:,ii,jj))
      CALL LeftEigenvectorsRotationalInvariancePrimitiveSystem(          W(1:nVar,ii,jj),NormVectY,LLY(:,:,ii,jj))
      CALL RightEigenvectorsRotationalInvariancePrimitiveSystem(         W(1:nVar,ii,jj),NormVectY,RRY(:,:,ii,jj))
      CALL TransitionMatrixConsToPrim(Cons(1:nVar),MCtoP(:,:,ii,jj))
      CALL TransitionMatrixPrimToCons(Cons(1:nVar),MPtoC(:,:,ii,jj))
  END DO
#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE ComputeTransitionMatricesPrimitiveInput
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ComputeTransitionMatricesConservedInput(W,vbXs,vbXe,vbYs,vbYe,lbXs,lbXe,lbYs,lbYe)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: LLX
USE MOD_FiniteVolume2D_vars,ONLY: RRX
USE MOD_FiniteVolume2D_vars,ONLY: LLY
USE MOD_FiniteVolume2D_vars,ONLY: RRY
USE MOD_FiniteVolume2D_vars,ONLY: MCtoP
USE MOD_FiniteVolume2D_vars,ONLY: MPtoC
USE MOD_FiniteVolume2D_vars,ONLY: NormVectX
USE MOD_FiniteVolume2D_vars,ONLY: NormVectY
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX, nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
#ifndef MULTIFLUID
USE MOD_Equation,           ONLY: LeftEigenvectorsXDirection
USE MOD_Equation,           ONLY: RightEigenvectorsXDirection
USE MOD_Equation,           ONLY: LeftEigenvectorsYDirection
USE MOD_Equation,           ONLY: RightEigenvectorsYDirection
#endif
USE MOD_Equation,           ONLY: TransitionMatrixConsToPrim
USE MOD_Equation,           ONLY: TransitionMatrixPrimToCons
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,vbXs:vbXe,vbYs:vbYe),    INTENT(IN)          :: W                   !*NB: Conserved variables assumed as inputs
INTEGER,                                        INTENT(IN)          :: vbXs,vbXe,vbYs,vbYe !*Bounds for the vector
INTEGER,                                        INTENT(IN)          :: lbXs,lbXe,lbYs,lbYe !*Bounds for the loop   
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, indi,indj
REAL    :: LMat(1:nVar,1:nVar)
REAL    :: RMat(1:nVar,1:nVar)
REAL    :: MMat(1:nVar,1:nVar)
!-------------------------------------------------------------------------------!


  DO ii=lbXs,lbXe
    DO jj=lbYs,lbYe
#ifndef MULTIFLUID
      CALL LeftEigenvectorsXDirection( W(1:nVar,ii,jj),LLX(:,:,ii,jj))
      CALL RightEigenvectorsXDirection(W(1:nVar,ii,jj),RRX(:,:,ii,jj))
      CALL LeftEigenvectorsYDirection( W(1:nVar,ii,jj),LLY(:,:,ii,jj))
      CALL RightEigenvectorsYDirection(W(1:nVar,ii,jj),RRY(:,:,ii,jj))
#endif
      CALL TransitionMatrixConsToPrim(W(1:nVar,ii,jj),MCtoP(:,:,ii,jj))
      CALL TransitionMatrixPrimToCons(W(1:nVar,ii,jj),MPtoC(:,:,ii,jj))
    END DO
  END DO 

  !*Ghost cells (only one layer)
  !*B,R,U,L

  !*B
  jj=lbYs-1
  DO ii=lbXs,lbXe
#ifndef MULTIFLUID
    CALL LeftEigenvectorsXDirection( W(1:nVar,ii,jj),LLX(:,:,ii,jj))
    CALL RightEigenvectorsXDirection(W(1:nVar,ii,jj),RRX(:,:,ii,jj))
    CALL LeftEigenvectorsYDirection( W(1:nVar,ii,jj),LLY(:,:,ii,jj))
    CALL RightEigenvectorsYDirection(W(1:nVar,ii,jj),RRY(:,:,ii,jj))
#endif
    CALL TransitionMatrixConsToPrim(W(1:nVar,ii,jj),MCtoP(:,:,ii,jj))
    CALL TransitionMatrixPrimToCons(W(1:nVar,ii,jj),MPtoC(:,:,ii,jj))
  END DO 

  !*R
  ii=lbXe+1
  DO jj=lbYs,lbYe
#ifndef MULTIFLUID
    CALL LeftEigenvectorsXDirection( W(1:nVar,ii,jj),LLX(:,:,ii,jj))
    CALL RightEigenvectorsXDirection(W(1:nVar,ii,jj),RRX(:,:,ii,jj))
    CALL LeftEigenvectorsYDirection( W(1:nVar,ii,jj),LLY(:,:,ii,jj))
    CALL RightEigenvectorsYDirection(W(1:nVar,ii,jj),RRY(:,:,ii,jj))
#endif
    CALL TransitionMatrixConsToPrim(W(1:nVar,ii,jj),MCtoP(:,:,ii,jj))
    CALL TransitionMatrixPrimToCons(W(1:nVar,ii,jj),MPtoC(:,:,ii,jj))
  END DO

  !*U
  jj=lbYe+1
  DO ii=lbXs,lbXe
#ifndef MULTIFLUID
    CALL LeftEigenvectorsXDirection( W(1:nVar,ii,jj),LLX(:,:,ii,jj))
    CALL RightEigenvectorsXDirection(W(1:nVar,ii,jj),RRX(:,:,ii,jj))
    CALL LeftEigenvectorsYDirection( W(1:nVar,ii,jj),LLY(:,:,ii,jj))
    CALL RightEigenvectorsYDirection(W(1:nVar,ii,jj),RRY(:,:,ii,jj))
#endif
    CALL TransitionMatrixConsToPrim(W(1:nVar,ii,jj),MCtoP(:,:,ii,jj))
    CALL TransitionMatrixPrimToCons(W(1:nVar,ii,jj),MPtoC(:,:,ii,jj))
  END DO 

  !*L
  ii=lbXs-1
  DO jj=lbYs,lbYe
#ifndef MULTIFLUID
    CALL LeftEigenvectorsXDirection( W(1:nVar,ii,jj),LLX(:,:,ii,jj))
    CALL RightEigenvectorsXDirection(W(1:nVar,ii,jj),RRX(:,:,ii,jj))
    CALL LeftEigenvectorsYDirection( W(1:nVar,ii,jj),LLY(:,:,ii,jj))
    CALL RightEigenvectorsYDirection(W(1:nVar,ii,jj),RRY(:,:,ii,jj))
#endif
    CALL TransitionMatrixConsToPrim(W(1:nVar,ii,jj),MCtoP(:,:,ii,jj))
    CALL TransitionMatrixPrimToCons(W(1:nVar,ii,jj),MPtoC(:,:,ii,jj))
  END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE ComputeTransitionMatricesConservedInput
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE FinalizeFiniteVolume()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: MeshNodes
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP  
USE MOD_FiniteVolume2D_vars,ONLY: NormVectX
USE MOD_FiniteVolume2D_vars,ONLY: NormVectY
USE MOD_FiniteVolume2D_vars,ONLY: TangVectX
USE MOD_FiniteVolume2D_vars,ONLY: TangVectY
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: Ut
USE MOD_FiniteVolume2D_vars,ONLY: UtWB
USE MOD_FiniteVolume2D_vars,ONLY: FYWB
USE MOD_FiniteVolume2D_vars,ONLY: FXWB
USE MOD_FiniteVolume2D_vars,ONLY: SWB
USE MOD_FiniteVolume2D_vars,ONLY: S
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: FX
USE MOD_FiniteVolume2D_vars,ONLY: FY
USE MOD_FiniteVolume2D_vars,ONLY: FluxX
USE MOD_FiniteVolume2D_vars,ONLY: FluxY
USE MOD_FiniteVolume2D_vars,ONLY: Ind

USE MOD_FiniteVolume2D_vars,ONLY: UN0
USE MOD_FiniteVolume2D_vars,ONLY: K0
USE MOD_FiniteVolume2D_vars,ONLY: K1
USE MOD_FiniteVolume2D_vars,ONLY: K2
USE MOD_FiniteVolume2D_vars,ONLY: K3
USE MOD_FiniteVolume2D_vars,ONLY: K4
USE MOD_FiniteVolume2D_vars,ONLY: K5
USE MOD_FiniteVolume2D_vars,ONLY: Ua
USE MOD_FiniteVolume2D_vars,ONLY: Up
USE MOD_FiniteVolume2D_vars,ONLY: FUp
USE MOD_FiniteVolume2D_vars,ONLY: Gravitational_Potential_Averages
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGP
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGPBnd

USE MOD_FiniteVolume2D_vars,ONLY: LLX
USE MOD_FiniteVolume2D_vars,ONLY: RRX
USE MOD_FiniteVolume2D_vars,ONLY: LLY
USE MOD_FiniteVolume2D_vars,ONLY: RRY
USE MOD_FiniteVolume2D_vars,ONLY: MCtoP
USE MOD_FiniteVolume2D_vars,ONLY: MPtoC

#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: WCt
USE MOD_FiniteVolume2D_vars,ONLY: UN0_WC
USE MOD_FiniteVolume2D_vars,ONLY: K0_WC
USE MOD_FiniteVolume2D_vars,ONLY: K1_WC
USE MOD_FiniteVolume2D_vars,ONLY: K2_WC
USE MOD_FiniteVolume2D_vars,ONLY: K3_WC
USE MOD_FiniteVolume2D_vars,ONLY: K4_WC
USE MOD_FiniteVolume2D_vars,ONLY: K5_WC
USE MOD_FiniteVolume2D_vars,ONLY: Wp_WC
USE MOD_FiniteVolume2D_vars,ONLY: Wa_WC
USE MOD_FiniteVolume2D_vars,ONLY: FWp_WC




!*The following structures are essential for IMEX, 
!*but I need them also for the explicit version, for the linear system to initialize the velocity in Gresho

USE MOD_FiniteVolume2D_vars,ONLY: From_2Indices_To_1Index_WC
USE MOD_FiniteVolume2D_vars,ONLY: From_1Index_To_2Indices_WC
USE MOD_FiniteVolume2D_vars,ONLY: Neighbours_1Index_WC


#ifdef IMEX
USE MOD_FiniteVolume2D_vars,ONLY: Values_WC
USE MOD_FiniteVolume2D_vars,ONLY: Columns_WC
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_WC
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_WC
USE MOD_FiniteVolume2D_vars,ONLY: rhs_WC
USE MOD_FiniteVolume2D_vars,ONLY: sol_WC
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_WC
#endif
#endif

#if defined(CENTEREDPRIMITIVE) && defined(PATHCONSERVATIVESHOCKDETECTION)
USE MOD_FiniteVolume2D_vars,ONLY: Cell_To_Limit_PCSD
#endif
#if defined(CENTEREDPRIMITIVE)
USE MOD_FiniteVolume2D_vars,ONLY: Errors_PCSD
#endif

#ifdef ACTIVEFLUX
!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary_X
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_X  
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: Wt_X

USE MOD_FiniteVolume2D_vars,ONLY: UN0_X
USE MOD_FiniteVolume2D_vars,ONLY: K0_X
USE MOD_FiniteVolume2D_vars,ONLY: K1_X
USE MOD_FiniteVolume2D_vars,ONLY: K2_X
USE MOD_FiniteVolume2D_vars,ONLY: K3_X
USE MOD_FiniteVolume2D_vars,ONLY: K4_X
USE MOD_FiniteVolume2D_vars,ONLY: K5_X
USE MOD_FiniteVolume2D_vars,ONLY: Wa_X
USE MOD_FiniteVolume2D_vars,ONLY: Wp_X
USE MOD_FiniteVolume2D_vars,ONLY: FWp_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary_Y
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_Y
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wt_Y

USE MOD_FiniteVolume2D_vars,ONLY: UN0_Y
USE MOD_FiniteVolume2D_vars,ONLY: K0_Y
USE MOD_FiniteVolume2D_vars,ONLY: K1_Y
USE MOD_FiniteVolume2D_vars,ONLY: K2_Y
USE MOD_FiniteVolume2D_vars,ONLY: K3_Y
USE MOD_FiniteVolume2D_vars,ONLY: K4_Y
USE MOD_FiniteVolume2D_vars,ONLY: K5_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wa_Y
USE MOD_FiniteVolume2D_vars,ONLY: Wp_Y
USE MOD_FiniteVolume2D_vars,ONLY: FWp_Y


USE MOD_FiniteVolume2D_vars,ONLY: W_qp_X
USE MOD_FiniteVolume2D_vars,ONLY: W_qp_Y

!*The following structures are essential for IMEX, 
!*but I need them also for the explicit version, for the linear system to initialize the velocity in Gresho
!*In any case, they are referred to W_X and W_Y and so they are meant to be present if ACTIVEFLUX flag is active

!*For staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: From_2Indices_To_1Index_X
USE MOD_FiniteVolume2D_vars,ONLY: From_1Index_To_2Indices_X
USE MOD_FiniteVolume2D_vars,ONLY: Neighbours_1Index_X

!*For staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: From_2Indices_To_1Index_Y
USE MOD_FiniteVolume2D_vars,ONLY: From_1Index_To_2Indices_Y
USE MOD_FiniteVolume2D_vars,ONLY: Neighbours_1Index_Y

#if defined(IMEX) || defined(IMEXMOMENTUM)
!*For staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_X
USE MOD_FiniteVolume2D_vars,ONLY: Columns_X
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_X
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_X
USE MOD_FiniteVolume2D_vars,ONLY: rhs_X
USE MOD_FiniteVolume2D_vars,ONLY: sol_X
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_X

!*For staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: Values_Y
USE MOD_FiniteVolume2D_vars,ONLY: Columns_Y
USE MOD_FiniteVolume2D_vars,ONLY: RowStart_Y
USE MOD_FiniteVolume2D_vars,ONLY: Diagonal_Y
USE MOD_FiniteVolume2D_vars,ONLY: rhs_Y
USE MOD_FiniteVolume2D_vars,ONLY: sol_Y
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations_Y
#endif

#ifdef MULTIFLUID
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_U
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_W_X
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Which_Fluid_In_Cell_U
#endif

#endif

#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
!*Shared
USE MOD_FiniteVolume2D_vars,ONLY: BX_Vol
USE MOD_FiniteVolume2D_vars,ONLY: BY_Vol
USE MOD_FiniteVolume2D_vars,ONLY: BX_Sur
USE MOD_FiniteVolume2D_vars,ONLY: BY_Sur
USE MOD_FiniteVolume2D_vars,ONLY: BX_Sur_qp
USE MOD_FiniteVolume2D_vars,ONLY: BY_Sur_qp
USE MOD_FiniteVolume2D_vars,ONLY: am_qp
USE MOD_FiniteVolume2D_vars,ONLY: ap_qp
#endif

#ifdef PATANKAR 
USE MOD_FiniteVolume2D_vars,ONLY: NNZsparse
USE MOD_FiniteVolume2D_vars,ONLY: ProductionSparse
USE MOD_FiniteVolume2D_vars,ONLY: DestructionSparse
USE MOD_FiniteVolume2D_vars,ONLY: RowStart
USE MOD_FiniteVolume2D_vars,ONLY: ColumnsVector
USE MOD_FiniteVolume2D_vars,ONLY: ProdUp
USE MOD_FiniteVolume2D_vars,ONLY: DestUp
USE MOD_FiniteVolume2D_vars,ONLY: SparseIndexMat
USE MOD_FiniteVolume2D_vars,ONLY: JacobiIterations
#endif   
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!

DEALLOCATE(MeshNodes)
DEALLOCATE(MeshBary)
DEALLOCATE(MeshGP)
DEALLOCATE(WeightsGP)
DEALLOCATE(WeightsGPBnd)


DEALLOCATE(U)
DEALLOCATE(Gravitational_Potential_Averages)
DEALLOCATE(V)
DEALLOCATE(Ut)
DEALLOCATE(WM)
DEALLOCATE(WP)
DEALLOCATE(S)
DEALLOCATE(FX)
DEALLOCATE(FY)
DEALLOCATE(FluxX)
DEALLOCATE(FluxY)
DEALLOCATE(UN0)
DEALLOCATE(K0)
DEALLOCATE(K1)
DEALLOCATE(K2)
DEALLOCATE(K3)
DEALLOCATE(K4)
DEALLOCATE(K5)
    
DEALLOCATE(Ua)
DEALLOCATE(Up)
DEALLOCATE(FUp)

DEALLOCATE(LLX) 
DEALLOCATE(RRX) 
DEALLOCATE(LLY) 
DEALLOCATE(RRY) 
DEALLOCATE(MCtoP) 
DEALLOCATE(MPtoC) 

#ifdef CENTEREDPRIMITIVE
DEALLOCATE( WC     )
DEALLOCATE( WCt    )
DEALLOCATE( UN0_WC )
DEALLOCATE( K0_WC  )
DEALLOCATE( K1_WC  )
DEALLOCATE( K2_WC  )
DEALLOCATE( K3_WC  )
DEALLOCATE( K4_WC  )
DEALLOCATE( K5_WC  )
DEALLOCATE( Wp_WC  )
DEALLOCATE( Wa_WC  )
DEALLOCATE( FWp_WC )

!*The following structures are essential for IMEX, 
!*but I need them also for the explicit version, for the linear system to initialize the velocity in Gresho

DEALLOCATE(From_2Indices_To_1Index_WC)
DEALLOCATE(From_1Index_To_2Indices_WC)
DEALLOCATE(Neighbours_1Index_WC)

#if defined(IMEX)
DEALLOCATE(Values_WC)
DEALLOCATE(Columns_WC)
DEALLOCATE(RowStart_WC)
DEALLOCATE(Diagonal_WC)
DEALLOCATE(rhs_WC)
DEALLOCATE(sol_WC)
#ifdef COUNTJACOBI
DEALLOCATE(JacobiIterations_WC)
#endif
#endif

#endif

#if defined(CENTEREDPRIMITIVE) && defined(PATHCONSERVATIVESHOCKDETECTION)
DEALLOCATE(Cell_To_Limit_PCSD)
#endif
#if defined(CENTEREDPRIMITIVE)
DEALLOCATE(Errors_PCSD)
#endif

#ifdef ACTIVEFLUX
!*Staggering in X direction
DEALLOCATE(MeshBary_X)
DEALLOCATE(MeshGP_X)
DEALLOCATE(W_X)
DEALLOCATE(Wt_X)

DEALLOCATE(UN0_X)
DEALLOCATE(K0_X)
DEALLOCATE(K1_X)
DEALLOCATE(K2_X)
DEALLOCATE(K3_X)
DEALLOCATE(K4_X)
DEALLOCATE(K5_X)
DEALLOCATE(Wa_X)
DEALLOCATE(Wp_X)
DEALLOCATE(FWp_X)

!*Staggering in y direction
DEALLOCATE(MeshBary_Y)
DEALLOCATE(MeshGP_Y)
DEALLOCATE(W_Y)
DEALLOCATE(Wt_Y)

DEALLOCATE(UN0_Y)
DEALLOCATE(K0_Y)
DEALLOCATE(K1_Y)
DEALLOCATE(K2_Y)
DEALLOCATE(K3_Y)
DEALLOCATE(K4_Y)
DEALLOCATE(K5_Y)
DEALLOCATE(Wa_Y)
DEALLOCATE(Wp_Y)
DEALLOCATE(FWp_Y)

DEALLOCATE(W_qp_X)
DEALLOCATE(W_qp_Y)

!*The following structures are essential for IMEX, 
!*but I need them also for the explicit version, for the linear system to initialize the velocity in Gresho
!*In any case, they are referred to W_X and W_Y and so they are meant to be present if ACTIVEFLUX flag is active

!*For staggering in X direction
DEALLOCATE(From_2Indices_To_1Index_X)
DEALLOCATE(From_1Index_To_2Indices_X)
DEALLOCATE(Neighbours_1Index_X)

!*For staggering in Y direction
DEALLOCATE(From_2Indices_To_1Index_Y)
DEALLOCATE(From_1Index_To_2Indices_Y)
DEALLOCATE(Neighbours_1Index_Y)


#if defined(IMEX) || defined(IMEXMOMENTUM)

!*For staggering in X direction
DEALLOCATE(Values_X)
DEALLOCATE(Columns_X)
DEALLOCATE(RowStart_X)
DEALLOCATE(Diagonal_X)
DEALLOCATE(rhs_X)
DEALLOCATE(sol_X)
#ifdef COUNTJACOBI
DEALLOCATE(JacobiIterations_X)
#endif

!*For staggering in Y direction
DEALLOCATE(Values_Y)
DEALLOCATE(Columns_Y)
DEALLOCATE(RowStart_Y)
DEALLOCATE(Diagonal_Y)
DEALLOCATE(rhs_Y)
DEALLOCATE(sol_Y)
#ifdef COUNTJACOBI
DEALLOCATE(JacobiIterations_Y)
#endif

#endif




#ifdef MULTIFLUID
DEALLOCATE(Troubled_Cell_U)
DEALLOCATE(Troubled_Cell_W_X)
DEALLOCATE(Troubled_Cell_W_Y)
DEALLOCATE(Which_Fluid_In_Cell_U)

#endif

#endif

#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
!*Shared
DEALLOCATE(BX_Vol)
DEALLOCATE(BY_Vol)
DEALLOCATE(BX_Sur)
DEALLOCATE(BY_Sur)
DEALLOCATE(BX_Sur_qp)
DEALLOCATE(BY_Sur_qp)
DEALLOCATE(am_qp)
DEALLOCATE(ap_qp)
#endif

#ifdef WELLBALANCED
DEALLOCATE(UtWB)
DEALLOCATE(FXWB)
DEALLOCATE(FYWB)
DEALLOCATE(SWB)
#endif

#ifdef PATANKAR
DEALLOCATE(ProductionSparse)
DEALLOCATE(DestructionSparse)
DEALLOCATE(ColumnsVector)
DEALLOCATE(RowStart)

DEALLOCATE(ProdUp)
DEALLOCATE(DestUp)
DEALLOCATE(SparseIndexMat)

DEALLOCATE(JacobiIterations)
#endif   

!-------------------------------------------------------------------------------!
END SUBROUTINE FinalizeFiniteVolume
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE FluxFX()
!-------------------------------------------------------------------------------!
USE MOD_Equation,           ONLY: NormalFlux1D
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: FX
USE MOD_FiniteVolume2D_vars,ONLY: FluxX
USE MOD_FiniteVolume2D_vars,ONLY: NormVectX, TangVectX
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGPBnd
USE MOD_FiniteVolume2D_vars,ONLY: W_qp_X
#ifdef MULTIFLUID
#ifdef RIEMANNPROBLEM
USE MOD_FiniteVolume2D_vars,ONLY: Which_Fluid_In_Cell_U
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE multifluid_riemann_mod, ONLY: multifluid_riemann
USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Equation,           ONLY: PrimToCons
USE MOD_Equation,           ONLY: Get_Gamma
USE MOD_Equation,           ONLY: Get_Pressure_Infinity
#endif
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, iGP
INTEGER :: indi, indj
#ifdef MULTIFLUID
#ifdef RIEMANNPROBLEM
REAL    :: UM(1:nVar), UP(1:nVar) !*U minus, U plus
REAL    :: VM(1:nVar), VP(1:nVar) !*V minus, V plus
INTEGER :: which_fluid
INTEGER :: count1,       count2
REAL    :: invdist1,     invdist2
REAL    :: UAvg1(1:nVar), UAvg2(1:nVar)
REAL    :: dist
INTEGER :: iCheck, jCheck
REAL    :: Gmm_L, pinf_L, density_L, velocity_L, pressure_L
REAL    :: Gmm_R, pinf_R, density_R, velocity_R, pressure_R
REAL    :: density_star_L, density_star_R, velocity_star, pressure_star
#endif
#endif
!-------------------------------------------------------------------------------!

FX    = 0.0
FluxX = 0.0

DO jj=1,nElemsY
  DO ii=0,nElemsX
    CALL NormalFlux1D(&
                      W_qp_X(1:nVar,1:nGPs,ii,jj),&
                      NormVectX(1:nDims),&
                      TangVectX(1:nDims),&
                      FluxX(1:nVar,1:nGPs,ii,jj))
  END DO
END DO

#ifdef MULTIFLUID
#ifdef RIEMANNPROBLEM
DO jj=1,nElemsY
  DO ii=1,nElemsX
    IF ( ABS(Which_Fluid_In_Cell_U(ii,jj)) .EQ. 2 ) THEN

      !*================================
      !*1) Replace left and right cell averages by single fluid ones
      !*================================

      !*Left
      iCheck=ii-1
      jCheck=jj
      which_fluid=Which_Fluid_In_Cell_U(iCheck,jCheck)

      IF (ABS(which_fluid) .EQ. 1) THEN !*The cell has a single fluid, we are done
        UM(1:nVar)=U(1:nVar,iCheck,jCheck)
      ELSE !*We need to replace the value

        !*Let us look for surrounding cells with single fluid
        count1       =0 
        invdist1     =0.0 
        UAvg1(1:nVar)=0.0 
        count2       =0 
        invdist2     =0.0 
        UAvg2(1:nVar)=0.0 


        DO indj=-1,1
          DO indi=-1,1

            !*Let us skip the cell itself
            IF ((indi .EQ. 0) .AND. (indj .EQ. 0)) THEN
              CYCLE
            END IF

            IF (ABS(indi)*ABS(indj).EQ.1) THEN !*Corners
              dist=SQRT(MESH_DX(1)**2+MESH_DX(2)**2)
            ELSE IF (ABS(indi) .EQ. 1) THEN !*right-left
              dist=MESH_DX(1)
            ELSE !*up-down
              dist=MESH_DX(2)
            END IF

            IF (Which_Fluid_In_Cell_U(iCheck+indi,jCheck+indj) .EQ. 1) THEN
              invdist1=invdist1+1.0/dist
              UAvg1=Uavg1+U(1:nVar,iCheck+indi,jCheck+indj)*1.0/dist
              count1=count1+1
            ELSE IF (Which_Fluid_In_Cell_U(iCheck+indi,jCheck+indj) .EQ. -1) THEN
              invdist2=invdist2+1.0/dist
              UAvg2=Uavg2+U(1:nVar,iCheck+indi,jCheck+indj)*1.0/dist
              count2=count2+1
            ELSE
              !*Do nothing
            END IF

          END DO
        END DO

        IF ((count1 .EQ. 0) .AND. (count2 .EQ. 0)) THEN
          PRINT*, "Problem, no single fluid cell found at left of a troubled cell"
          STOP
        ELSE
          IF (which_fluid .EQ. 2) THEN
            IF (count1.GT.0) THEN
              UM=Uavg1*invdist1
              UM(nVar)=1.0
            ELSE
              UM=Uavg2*invdist2
              UM(nVar)=-1.0
            END IF
          ELSE
            IF (count2.GT.0) THEN
              UM=Uavg2*invdist2
              UM(nVar)=-1.0
            ELSE
              UM=Uavg1*invdist1
              UM(nVar)=1.0
            END IF
          END IF

        END IF

      END IF

      !*-----------------------------------
      !*PRINT*, UM-U(1:nVar, iCheck,jCheck)
      !*-----------------------------------

      !*Right
      iCheck=ii+1
      jCheck=jj
      which_fluid=Which_Fluid_In_Cell_U(iCheck,jCheck)

      IF (ABS(which_fluid) .EQ. 1) THEN !*The cell has a single fluid, we are done
        UP(1:nVar)=U(1:nVar,iCheck,jCheck)
      ELSE !*We need to replace the value

        !*Let us look for surrounding cells with single fluid
        count1       =0 
        invdist1     =0.0 
        UAvg1(1:nVar)=0.0 
        count2       =0 
        invdist2     =0.0 
        UAvg2(1:nVar)=0.0 


        DO indj=-1,1
          DO indi=-1,1

            !*Let us skip the cell itself
            IF ((indi .EQ. 0) .AND. (indj .EQ. 0)) THEN
              CYCLE
            END IF

            IF (ABS(indi)*ABS(indj).EQ.1) THEN !*Corners
              dist=SQRT(MESH_DX(1)**2+MESH_DX(2)**2)
            ELSE IF (ABS(indi) .EQ. 1) THEN !*right-left
              dist=MESH_DX(1)
            ELSE !*up-down
              dist=MESH_DX(2)
            END IF

            IF (Which_Fluid_In_Cell_U(iCheck+indi,jCheck+indj) .EQ. 1) THEN
              invdist1=invdist1+1.0/dist
              UAvg1=Uavg1+U(1:nVar,iCheck+indi,jCheck+indj)*1.0/dist
              count1=count1+1
            ELSE IF (Which_Fluid_In_Cell_U(iCheck+indi,jCheck+indj) .EQ. -1) THEN
              invdist2=invdist2+1.0/dist
              UAvg2=Uavg2+U(1:nVar,iCheck+indi,jCheck+indj)*1.0/dist
              count2=count2+1
            ELSE
              !*Do nothing
            END IF

          END DO
        END DO

        IF ((count1 .EQ. 0) .AND. (count2 .EQ. 0)) THEN
          PRINT*, "Problem, no single fluid cell found at right of a troubled cell"
          STOP
        ELSE
          IF (which_fluid .EQ. 2) THEN
            IF (count1.GT.0) THEN
              UP=Uavg1*invdist1
              UP(nVar)=1.0
            ELSE
              UP=Uavg2*invdist2
              UP(nVar)=-1.0
            END IF
          ELSE
            IF (count2.GT.0) THEN
              UP=Uavg2*invdist2
              UP(nVar)=-1.0
            ELSE
              UP=Uavg1*invdist1
              UP(nVar)=1.0
            END IF
          END IF

        END IF

      END IF

      !*-----------------------------------
      ! *PRINT*, UP-U(1:nVar, iCheck,jCheck)
      !*-----------------------------------



      !*================================
      !*2) Replace fluxes using the Multilfuid RP
      !*================================
      CALL ConsToPrim(UM(1:nVar),VM(1:nVar))
      CALL ConsToPrim(UP(1:nVar),VP(1:nVar))

      density_L  = VM(1)
      velocity_L = VM(2)
      pressure_L = VM(4)
      Gmm_L =Get_Gamma(VM(nVar))
      pinf_L=Get_Pressure_Infinity(VM(nVar))      

      density_R  = VP(1)
      velocity_R = VP(2)
      pressure_R = VP(4)
      Gmm_R =Get_Gamma(VP(nVar))
      pinf_R=Get_Pressure_Infinity(VP(nVar))      

      density_star_L = 0.0
      density_star_R = 0.0
      velocity_star  = 0.0
      pressure_star  = 0.0
      ! CALL multifluid_riemann(Gmm_L, pinf_L, density_L, velocity_L, pressure_L, Gmm_R, pinf_R, density_R, velocity_R, pressure_R, density_star_L, density_star_R, velocity_star, pressure_star)










    END IF
  END DO
END DO
#endif
#endif

DO iGP = 1,nGPs
  FX(:,:,:) = FX(:,:,:) + weightsGPBnd(iGP)*FluxX(:,iGP,:,:)
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE FluxFX
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE FluxFY()
!-------------------------------------------------------------------------------!
USE MOD_Equation,           ONLY: NormalFlux1D
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: FY
USE MOD_FiniteVolume2D_vars,ONLY: FluxY
USE MOD_FiniteVolume2D_vars,ONLY: NormVectY, TangVectY
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGPBnd
USE MOD_FiniteVolume2D_vars,ONLY: W_qp_Y
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, iGP
!-------------------------------------------------------------------------------!

FY    = 0.0
FluxY = 0.0

DO jj=0,nElemsY
  DO ii=1,nElemsX
    CALL NormalFlux1D(&
                       W_qp_Y(1:nVar,1:nGPs,ii,jj),&
                       NormVectY(1:nDims),&
                       TangVectY(1:nDims),&
                       FluxY(1:nVar,1:nGPs,ii,jj))
  END DO
END DO

DO iGP = 1,nGPs
  FY(:,:,:) = FY(:,:,:) + weightsGPBnd(iGP)*FluxY(:,iGP,:,:)
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE FluxFY
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE Adjust_Interface_values()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: qq, ii, jj, iGP
!-------------------------------------------------------------------------------!


#ifdef MOMENTUMINPRIMITIVEVARIABLES
#ifdef RECONSTRUCTIONVELOCITY
DO qq=1,nGPs
  DO jj=0,nElemsY+1+1
    DO ii=0,nElemsX+1+1
      IF (WM(1,qq,ii,jj) .NE. 0.0) THEN
        WM(2:3,qq,ii,jj)=WM(2:3,qq,ii,jj)*WM(1,qq,ii,jj)
      END IF
      IF (WP(1,qq,ii,jj) .NE. 0.0) THEN
        WP(2:3,qq,ii,jj)=WP(2:3,qq,ii,jj)*WP(1,qq,ii,jj)
      END IF
    END DO
  END DO
END DO
#endif
#endif

#ifdef MULTIFLUID
#ifdef RECONSTRUCTIONPHI
DO qq=1,nGPs
  DO jj=0,nElemsY+1+1
    DO ii=0,nElemsX+1+1
      IF (WM(1,qq,ii,jj) .NE. 0.0) THEN
        WM(nVar,qq,ii,jj)=WM(nVar,qq,ii,jj)*WM(1,qq,ii,jj)
      END IF
      IF (WP(1,qq,ii,jj) .NE. 0.0) THEN
        WP(nVar,qq,ii,jj)=WP(nVar,qq,ii,jj)*WP(1,qq,ii,jj)
      END IF
    END DO
  END DO
END DO
#endif
#endif


!-------------------------------------------------------------------------------!
END SUBROUTINE Adjust_Interface_values
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE NumericalFluxFX_W_INPUT(lbXs,lbXe,lbYs,lbYe)
!-------------------------------------------------------------------------------!
USE MOD_Equation,           ONLY: NumericalFluxPrimitiveSystem
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: FX
USE MOD_FiniteVolume2D_vars,ONLY: FluxX
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: NormVectX, TangVectX
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGPBnd
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
INTEGER,  INTENT(IN) :: lbXs,lbXe,lbYs,lbYe !*Bounds for the loop
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, iGP
!-------------------------------------------------------------------------------!

FX    = 0.0
FluxX = 0.0

DO jj=lbYs,lbYe
  DO ii=lbXs,lbXe
#ifdef MULTIFLUID
    CALL NumericalFluxPrimitiveSystem(&
                       WP(1:nVar+1,1:nGPs,ii+0,jj),&
                       WM(1:nVar+1,1:nGPs,ii+1,jj),&
                       NormVectX(1:nDims),&
                       TangVectX(1:nDims),&
                       FluxX(1:nVar,1:nGPs,ii,jj))
#else
    CALL NumericalFluxPrimitiveSystem(&
                       WP(1:nVar,1:nGPs,ii+0,jj),&
                       WM(1:nVar,1:nGPs,ii+1,jj),&
                       NormVectX(1:nDims),&
                       TangVectX(1:nDims),&
                       FluxX(1:nVar,1:nGPs,ii,jj))
#endif
  END DO
END DO

DO iGP = 1,nGPs
  FX(:,:,:) = FX(:,:,:) + weightsGPBnd(iGP)*FluxX(:,iGP,:,:)
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE NumericalFluxFX_W_INPUT
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE NumericalFluxFY_W_INPUT(lbXs,lbXe,lbYs,lbYe)
!-------------------------------------------------------------------------------!
USE MOD_Equation,           ONLY: NumericalFluxPrimitiveSystem
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: FY
USE MOD_FiniteVolume2D_vars,ONLY: FluxY
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: NormVectY, TangVectY
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGPBnd
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
INTEGER,  INTENT(IN) :: lbXs,lbXe,lbYs,lbYe !*Bounds for the loop
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, iGP
!-------------------------------------------------------------------------------!

FY    = 0.0
FluxY = 0.0

DO jj=lbYs,lbYe
  DO ii=lbXs,lbXe
#ifdef MULTIFLUID
    CALL NumericalFluxPrimitiveSystem(&
                       WP(1:nVar+1,1:nGPs,ii,jj+0),&
                       WM(1:nVar+1,1:nGPs,ii,jj+1),&
                       NormVectY(1:nDims),&
                       TangVectY(1:nDims),&
                       FluxY(1:nVar,1:nGPs,ii,jj))
#else
    CALL NumericalFluxPrimitiveSystem(&
                       WP(1:nVar,1:nGPs,ii,jj+0),&
                       WM(1:nVar,1:nGPs,ii,jj+1),&
                       NormVectY(1:nDims),&
                       TangVectY(1:nDims),&
                       FluxY(1:nVar,1:nGPs,ii,jj))
#endif
  END DO
END DO

DO iGP = 1,nGPs
  FY(:,:,:) = FY(:,:,:) + weightsGPBnd(iGP)*FluxY(:,iGP,:,:)
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE NumericalFluxFY_W_INPUT
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE BSurfX_W_INPUT(lbXs,lbXe,lbYs,lbYe)
!-------------------------------------------------------------------------------!
USE MOD_Equation,           ONLY: PathConservativeSurfaceContribution
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: BX_Sur    !*Equivalent for FX (with a dimension more for the side of the edge to take into account the local speed of propagation)
USE MOD_FiniteVolume2D_vars,ONLY: BX_Sur_qp !*Equivalent for FluxX
USE MOD_FiniteVolume2D_vars,ONLY: am_qp    
USE MOD_FiniteVolume2D_vars,ONLY: ap_qp    
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: NormVectX, TangVectX
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGPBnd
USE MOD_FiniteVolume2D_vars,ONLY: EPS_AVOIDING_DIVISION_BY_ZERO_IN_CU
USE MOD_FiniteVolume2D_vars,ONLY: WhichRiemannSolverPrimitiveSystem
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
INTEGER,  INTENT(IN) :: lbXs,lbXe,lbYs,lbYe !*Bounds for the loop
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, iGP
!-------------------------------------------------------------------------------!

BX_Sur    = 0.0
BX_Sur_qp = 0.0
am_qp     = 0.0
ap_qp     = 0.0

DO jj=lbYs,lbYe
  DO ii=lbXs,lbXe
#ifdef MULTIFLUID
    CALL PathConservativeSurfaceContribution(&
                       WP(1:nVar+1,1:nGPs,ii+0,jj),&
                       WM(1:nVar+1,1:nGPs,ii+1,jj),&
                       NormVectX(1:nDims),&
                       TangVectX(1:nDims),&
                       BX_Sur_qp(1:nVar,1:nGPs,ii,jj),&
                       am_qp(1:nGPs,ii,jj),&
                       ap_qp(1:nGPs,ii,jj))
#else
    CALL PathConservativeSurfaceContribution(&
                       WP(1:nVar,1:nGPs,ii+0,jj),&
                       WM(1:nVar,1:nGPs,ii+1,jj),&
                       NormVectX(1:nDims),&
                       TangVectX(1:nDims),&
                       BX_Sur_qp(1:nVar,1:nGPs,ii,jj),&
                       am_qp(1:nGPs,ii,jj),&
                       ap_qp(1:nGPs,ii,jj))
#endif
  END DO
END DO

!*NB: Here I cannot simply integrate. I have to keep into account ap and am

DO jj=lbYs,lbYe
  DO ii=lbXs,lbXe
    DO iGP = 1,nGPs
      IF ( WhichRiemannSolverPrimitiveSystem .EQ. 0 ) THEN !*Rusanov
        BX_Sur(:,ii,jj,1) = BX_Sur(:,ii,jj,1) + weightsGPBnd(iGP)*BX_Sur_qp(:,iGP,ii,jj)*(-0.5) !*Left  1 (to be used as right contribution by ii  )
        BX_Sur(:,ii,jj,2) = BX_Sur(:,ii,jj,2) + weightsGPBnd(iGP)*BX_Sur_qp(:,iGP,ii,jj)*( 0.5) !*Right 2 (to be used as left  contribution by ii+1)
      ELSE
        BX_Sur(:,ii,jj,1) = BX_Sur(:,ii,jj,1) + weightsGPBnd(iGP)*BX_Sur_qp(:,iGP,ii,jj)*am_qp(iGP,ii,jj)/(ap_qp(iGP,ii,jj)-am_qp(iGP,ii,jj)) !*Left  1 (to be used as right contribution by ii  )
        BX_Sur(:,ii,jj,2) = BX_Sur(:,ii,jj,2) + weightsGPBnd(iGP)*BX_Sur_qp(:,iGP,ii,jj)*ap_qp(iGP,ii,jj)/(ap_qp(iGP,ii,jj)-am_qp(iGP,ii,jj)) !*Right 2 (to be used as left  contribution by ii+1)
      END IF
    END DO
  END DO
END DO

!*Reference
!*PATH-CONSERVATIVE CENTRAL-UPWIND SCHEMES FOR NONCONSERVATIVE HYPERBOLIC SYSTEMS, Manuel Jesús Castro Dı́az, Alexander Kurganov, Tomás Morales de Luna
!*Eq (3.8)

!-------------------------------------------------------------------------------!
END SUBROUTINE BSurfX_W_INPUT
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE BSurfY_W_INPUT(lbXs,lbXe,lbYs,lbYe)
!-------------------------------------------------------------------------------!
USE MOD_Equation,           ONLY: PathConservativeSurfaceContribution
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: BY_Sur    !*Equivalent for FY (with a dimension more for the side of the edge to take into account the local speed of propagation)
USE MOD_FiniteVolume2D_vars,ONLY: BY_Sur_qp !*Equivalent for FluxY
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: NormVectY, TangVectY
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGPBnd
USE MOD_FiniteVolume2D_vars,ONLY: am_qp    
USE MOD_FiniteVolume2D_vars,ONLY: ap_qp    
USE MOD_FiniteVolume2D_vars,ONLY: EPS_AVOIDING_DIVISION_BY_ZERO_IN_CU
USE MOD_FiniteVolume2D_vars,ONLY: WhichRiemannSolverPrimitiveSystem
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
INTEGER,  INTENT(IN) :: lbXs,lbXe,lbYs,lbYe !*Bounds for the loop
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, iGP
!-------------------------------------------------------------------------------!

BY_Sur    = 0.0
BY_Sur_qp = 0.0
am_qp     = 0.0
ap_qp     = 0.0

DO jj=lbYs,lbYe
  DO ii=lbXs,lbXe
#ifdef MULTIFLUID
    CALL PathConservativeSurfaceContribution(&
                       WP(1:nVar+1,1:nGPs,ii,jj+0),&
                       WM(1:nVar+1,1:nGPs,ii,jj+1),&
                       NormVectY(1:nDims),&
                       TangVectY(1:nDims),&
                       BY_Sur_qp(1:nVar,1:nGPs,ii,jj),&
                       am_qp(1:nGPs,ii,jj),&
                       ap_qp(1:nGPs,ii,jj))
#else
    CALL PathConservativeSurfaceContribution(&
                       WP(1:nVar,1:nGPs,ii,jj+0),&
                       WM(1:nVar,1:nGPs,ii,jj+1),&
                       NormVectY(1:nDims),&
                       TangVectY(1:nDims),&
                       BY_Sur_qp(1:nVar,1:nGPs,ii,jj),&
                       am_qp(1:nGPs,ii,jj),&
                       ap_qp(1:nGPs,ii,jj))
#endif
  END DO
END DO

!*NB: Here I cannot simply integrate. I have to keep into account ap and am

DO jj=lbYs,lbYe
  DO ii=lbXs,lbXe
    DO iGP = 1,nGPs
      IF ( WhichRiemannSolverPrimitiveSystem .EQ. 0 ) THEN !*Rusanov
        BY_Sur(:,ii,jj,1) = BY_Sur(:,ii,jj,1) + weightsGPBnd(iGP)*BY_Sur_qp(:,iGP,ii,jj)*(-0.5) !*Down  1 (to be used as up    contribution by jj  )
        BY_Sur(:,ii,jj,2) = BY_Sur(:,ii,jj,2) + weightsGPBnd(iGP)*BY_Sur_qp(:,iGP,ii,jj)*( 0.5) !*Up    2 (to be used as down  contribution by jj+1)
      ELSE
        BY_Sur(:,ii,jj,1) = BY_Sur(:,ii,jj,1) + weightsGPBnd(iGP)*BY_Sur_qp(:,iGP,ii,jj)*am_qp(iGP,ii,jj)/(ap_qp(iGP,ii,jj)-am_qp(iGP,ii,jj)) !*Down  1 (to be used as up    contribution by jj  )
        BY_Sur(:,ii,jj,2) = BY_Sur(:,ii,jj,2) + weightsGPBnd(iGP)*BY_Sur_qp(:,iGP,ii,jj)*ap_qp(iGP,ii,jj)/(ap_qp(iGP,ii,jj)-am_qp(iGP,ii,jj)) !*Up    2 (to be used as down  contribution by jj+1)
      END IF
    END DO
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE BSurfY_W_INPUT
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE BVolX_W_INPUT(W,vbXs,vbXe,vbYs,vbYe,lbXs,lbXe,lbYs,lbYe,ReconstructionInput,Troubled_Cell_W)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: BX_Vol
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGP
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_Reconstruction     ,ONLY: MUSCL
USE MOD_Reconstruction     ,ONLY: kMUSCL
USE MOD_Reconstruction     ,ONLY: COMUSCL
USE MOD_Reconstruction     ,ONLY: VLMUSCL
USE MOD_Reconstruction     ,ONLY: M2MUSCL
USE MOD_Reconstruction     ,ONLY: VAMUSCL
USE MOD_Reconstruction     ,ONLY: SBMMUSCL
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
#ifdef IMEX
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: min_p
#endif
USE MOD_FiniteVolume2D_vars,ONLY: timescheme
#ifdef MULTIFLUID
USE MOD_Equation           ,ONLY: Get_Gamma
USE MOD_Equation           ,ONLY: Get_Pressure_Infinity
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,vbXs:vbXe,vbYs:vbYe), INTENT(IN)           :: W
INTEGER,                                     INTENT(IN)           :: vbXs,vbXe,vbYs,vbYe !*Bounds for the vector
INTEGER,                                     INTENT(IN)           :: lbXs,lbXe,lbYs,lbYe !*Bounds for the loop                                
INTEGER, DIMENSION(vbXs:vbXe,vbYs:vbYe),     INTENT(IN), OPTIONAL :: Troubled_Cell_W
INTEGER,                                     INTENT(IN), OPTIONAL :: ReconstructionInput
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER          :: ReconstructionUsed
REAL             :: BVol_in_qp(1:nVar,nGPs,nGPs,nElemsX+1,nElemsY+1)  !*NB: +1 in X direction and also in Y to handle both W_X and W_Y
REAL             :: VtempM(1:nVar,1:nGPs)
REAL             :: VtempP(1:nVar,1:nGPs)
REAL             :: VtempMsupport(1:nVar,1:nGPs)
REAL             :: VtempPsupport(1:nVar,1:nGPs)
REAL             :: dxp, dxvy, dxvx, dxro, dxrovx
INTEGER          :: ii, jj, iVar, iGP, jGP
CHARACTER(LEN=255) :: ErrorMessage
#ifdef MULTIFLUID
REAL             :: roPhi
REAL             :: dxPhi
REAL             :: pinf
#endif
REAL             :: ro, vx, vy, p
!-------------------------------------------------------------------------------!

IF (PRESENT(ReconstructionInput)) THEN
  ReconstructionUsed=ReconstructionInput
ELSE
  ReconstructionUsed=Reconstruction
END IF

BX_Vol     = 0.0 
BVol_in_qp = 0.0


DO jj=lbYs,lbYe
  DO ii=lbXs,lbXe

#ifdef PWCINTROUBLEDCELLS
#ifdef MULTIFLUID
    IF(PRESENT(Troubled_Cell_W)) THEN
      IF (Troubled_Cell_W(ii,jj) .NE. 0) THEN !*PWC
        CYCLE
      END IF
    END IF
#else
PRINT*, "Impossible to use PWCINTROUBLEDCELLS"
PRINT*, "Disable this flag"
PRINT*, "You should not be here"
STOP
#endif
#endif

    SELECT CASE (ABS(ReconstructionUsed))
      CASE(1)
        !*Do nothing. Keep BX_Vol = 0.0 
        !*Because the reconstruction is PWC and the derivatives disappear
        !*NB: Consistency is ensured by the surface contribution of PC
    
      CASE(10,2,20,21,22,23,24,25,26,27,28,29) !*MUSCL, k-MUSCL, MUSCL_CO, MUSCL_VL, MUSCL_M, MUSCL_VA, MUSCL_SBM

        VtempMsupport(:,nGPs)=WM(1:nVar,nGPs,ii,jj)
        VtempPsupport(:,nGPs)=WP(1:nVar,nGPs,ii,jj)

#if MOMENTUMINPRIMITIVEVARIABLES
        VtempMsupport(2:3,nGPs)=VtempMsupport(2:3,nGPs)/VtempMsupport(1,nGPs)
        VtempPsupport(2:3,nGPs)=VtempPsupport(2:3,nGPs)/VtempPsupport(1,nGPs)
#endif

        VtempM(1:nVar,nGPs)=VtempMsupport(1:nVar,nGPs)
        VtempP(1:nVar,nGPs)=VtempPsupport(1:nVar,nGPs)

        !*ro 1; vx 2; vy 3; p  4
        dxro = (VtempP(1,nGPs)-VtempM(1,nGPs))

        dxvx = (VtempP(2,nGPs)-VtempM(2,nGPs))
        dxvy = (VtempP(3,nGPs)-VtempM(3,nGPs))

        dxp  = (VtempP(4,nGPs)-VtempM(4,nGPs))
#ifdef MULTIFLUID
#ifdef PRIMITIVEFORMULATIONLEVELSET
        dxPhi= (VtempP(5,nGPs)-VtempM(5,nGPs))
#endif
#endif

        dxrovx = (VtempP(1,nGPs)*VtempP(2,nGPs)-VtempM(1,nGPs)*VtempM(2,nGPs))


        ro = W(1,ii,jj)
        p  = W(4,ii,jj)
#if defined(MOMENTUMINPRIMITIVEVARIABLES)
        !*In this case I am entering with momentum
        vx = W(2,ii,jj)/W(1,ii,jj)
        vy = W(3,ii,jj)/W(1,ii,jj)
#else
        vx = W(2,ii,jj)
        vy = W(3,ii,jj)
#endif


        BVol_in_qp(1,nGPs,nGPs,ii,jj) =  0.

#ifndef MOMENTUMINPRIMITIVEVARIABLES
!*START not MOMENTUMINPRIMITIVEVARIABLES

#ifdef IMEX
  !*START IMEX

#ifdef IMEXL2FULLYUPWINDED
    !*START IMEXL2FULLYUPWINDED
    BVol_in_qp(2,nGPs,nGPs,ii,jj) = -     1.0/ro * dxp
#else
    !*ELSE IMEXL2FULLYUPWINDED
    BVol_in_qp(2,nGPs,nGPs,ii,jj) = -     (max_ro-ro)/(ro*max_ro) * dxp
#endif
    !*END IMEXL2FULLYUPWINDED

#else
  !*ELSE IMEX

#if(1==1)
  BVol_in_qp(2,nGPs,nGPs,ii,jj) = -     1.0/ro * dxp                          !*What I suggest to use 
#else
  BVol_in_qp(2,nGPs,nGPs,ii,jj) = -0.5*(1.0/VtempP(1,nGPs)+1.0/VtempM(1,nGPs)) * dxp  !*What Alex initially proposed
#endif

#endif
  !*END IMEX

#ifdef RELAXATION
  !*START RELAXATION
  BVol_in_qp(2,nGPs,nGPs,ii,jj) = BVol_in_qp(2,nGPs,nGPs,ii,jj)/EPS_LM**2 !*In case I divide by EPS_LM**2
#endif
  !*END RELAXATION
  BVol_in_qp(3,nGPs,nGPs,ii,jj) = -         vx * dxvy   !*NB: W(2,ii,jj) is already the average
#endif
!*END not MOMENTUMINPRIMITIVEVARIABLES

        !*=======================
        !*SPLITTING
        !*=======================
#ifdef MULTIFLUID
roPhi=W(nVar,ii,jj)
Gmm=Get_Gamma(roPhi)
pinf=Get_Pressure_Infinity(roPhi)
#endif


#if defined(ALTERNATIVEFORMULATIONPRESSURE)
!*START ALTERNATIVEFORMULATIONPRESSURE
! PRINT*, "This is a consistent formulation."
! PRINT*, "It works but it is not suitable for the splitting."
! STOP

#ifdef MULTIFLUID
  !*START MULTIFLUID
  BVol_in_qp(4,nGPs,nGPs,ii,jj) = -( (Gmm-1.0)*p + Gmm*pinf ) *dxvx
#else
  !*ELSE MULTIFLUID
  BVol_in_qp(4,nGPs,nGPs,ii,jj) = -( (Gmm-1.0)*p ) *dxvx
#endif
  !*END MULTIFLUID 

#ifdef IMEX
  !*START IMEX
  PRINT*, "In particular, for IMEX one really can't do the splitting in this case."
  STOP
#endif
  !*END IMEX

#elif defined(PRESSUREFORMULATIONFORIMEXMOMENTUM)
!*ELSE ALTERNATIVEFORMULATIONPRESSURE

#ifdef IMEXMOMENTUM
  !*START IMEXMOMENTUM
  BVol_in_qp(4,nGPs,nGPs,ii,jj) = Gmm*p/ro*vx * dxro &
                                & -vx *dxp
#else
  !*ELSE IMEXMOMENTUM
  BVol_in_qp(4,nGPs,nGPs,ii,jj) = Gmm*p/ro*vx * dxro &
                                & -Gmm*p/vx* dxrovx &
                                & -vx *dxp
#endif
  !*END IMEXMOMENTUM

#else
!*ELSE ALTERNATIVEFORMULATIONPRESSURE

#ifdef IMEX
  !*START IMEX

#ifdef IMEXL2FULLYUPWINDED
    !*START IMEXL2FULLYUPWINDED
    BVol_in_qp(4,nGPs,nGPs,ii,jj) = - Gmm * p * dxvx &
                                      & -       vx * dxp
#else
    !*ELSE IMEXL2FULLYUPWINDED
    BVol_in_qp(4,nGPs,nGPs,ii,jj) = - Gmm * (p-min_p) * dxvx & !*This must be splitted 
                                      & -       vx * dxp            !*NB: This other one not
#endif
    !*END IMEXL2FULLYUPWINDED

#else
  !*ELSE IMEX


#ifdef MULTIFLUID
    !*START MULTIFLUID
    BVol_in_qp(4,nGPs,nGPs,ii,jj) = - Gmm * (p+pinf) * dxvx &
                                  & -       vx * dxp
#else
    !*ELSE MULTIFLUID
    BVol_in_qp(4,nGPs,nGPs,ii,jj) = - Gmm * p * dxvx &
                                  & -       vx * dxp
#endif
    !*END MULTIFLUID

#ifdef MULTIFLUID
    !*START MULTIFLUID

#ifndef PRIMITIVEFORMULATIONLEVELSET
      !*START PRIMITIVEFORMULATIONLEVELSET
      BVol_in_qp(nVar,nGPs,nGPs,ii,jj) = 0.0
#else
      !*ELSE PRIMITIVEFORMULATIONLEVELSET
      BVol_in_qp(nVar,nGPs,nGPs,ii,jj) = -vx * dxPhi
#endif
      !*END PRIMITIVEFORMULATIONLEVELSET

#endif
    !*END MULTIFLUID

#endif
  !*END IMEX

#endif
!*END ALTERNATIVEFORMULATIONPRESSURE


        !*--------------------------
        !*SAFETY CHECK
        !*--------------------------
        ! IF (ABS(W(2,ii,jj)-0.5*(VtempP(2,nGPs)+VtempM(2,nGPs))) .GT. 1e-15) THEN
        !   PRINT*, "Problem in PC volume X contribution"
        !   PRINT*, W(2,ii,jj)-0.5*(VtempP(2,nGPs)+VtempM(2,nGPs))
        !   STOP
        ! END IF
        !*--------------------------

      CASE DEFAULT
        ErrorMessage = "Reconstruction not implemented in BVolX_W_INPUT"
        WRITE(*,*) ErrorMessage
        STOP
    END SELECT
  END DO
END DO


DO jj=lbYs,lbYe
  DO ii=lbXs,lbXe
    DO iGP=1,nGPs
      DO jGP=1,nGPs
        BX_Vol(1:nVar,ii,jj) = BX_Vol(1:nVar,ii,jj) + WeightsGP(iGP,jGP) * BVol_in_qp(1:nVar,iGP,jGP,ii,jj) 
      END DO 
    END DO 
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE BVolX_W_INPUT
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE BVolY_W_INPUT(W,vbXs,vbXe,vbYs,vbYe,lbXs,lbXe,lbYs,lbYe,ReconstructionInput,Troubled_Cell_W)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: BY_Vol
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
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_Reconstruction     ,ONLY: MUSCL
USE MOD_Reconstruction     ,ONLY: kMUSCL
USE MOD_Reconstruction     ,ONLY: COMUSCL
USE MOD_Reconstruction     ,ONLY: VLMUSCL
USE MOD_Reconstruction     ,ONLY: M2MUSCL
USE MOD_Reconstruction     ,ONLY: VAMUSCL
USE MOD_Reconstruction     ,ONLY: SBMMUSCL
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
#ifdef IMEX
USE MOD_FiniteVolume2D_vars,ONLY: max_ro
USE MOD_FiniteVolume2D_vars,ONLY: min_p
#endif
USE MOD_FiniteVolume2D_vars,ONLY: timescheme
#ifdef MULTIFLUID
USE MOD_Equation           ,ONLY: Get_Gamma
USE MOD_Equation           ,ONLY: Get_Pressure_Infinity
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,vbXs:vbXe,vbYs:vbYe), INTENT(IN)           :: W
INTEGER,                                     INTENT(IN)           :: vbXs,vbXe,vbYs,vbYe !*Bounds for the vector
INTEGER,                                     INTENT(IN)           :: lbXs,lbXe,lbYs,lbYe !*Bounds for the loop   
INTEGER, DIMENSION(vbXs:vbXe,vbYs:vbYe),     INTENT(IN), OPTIONAL :: Troubled_Cell_W
INTEGER,                                     INTENT(IN), OPTIONAL :: ReconstructionInput
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER          :: ReconstructionUsed
REAL             :: BVol_in_qp(1:nVar,nGPs,nGPs,nElemsX+1,nElemsY+1)  !*NB: +1 in X direction and also in Y to handle both W_X and W_Y
REAL             :: VtempM(1:nVar,1:nGPs)
REAL             :: VtempP(1:nVar,1:nGPs)
REAL             :: VtempMsupport(1:nVar,1:nGPs)
REAL             :: VtempPsupport(1:nVar,1:nGPs)
REAL             :: dyp, dyvx, dyvy, dyro, dyrovy
INTEGER          :: ii, jj, iVar, iGP, jGP
CHARACTER(LEN=255) :: ErrorMessage
#ifdef MULTIFLUID
REAL             :: roPhi
REAL             :: dyPhi
REAL             :: pinf
#endif
REAL             :: ro, vx, vy, p
!-------------------------------------------------------------------------------!

IF (PRESENT(ReconstructionInput)) THEN
  ReconstructionUsed=ReconstructionInput
ELSE
  ReconstructionUsed=Reconstruction
END IF

BY_Vol     = 0.0 
BVol_in_qp = 0.0

DO jj=lbYs,lbYe
  DO ii=lbXs,lbXe

#ifdef PWCINTROUBLEDCELLS
#ifdef MULTIFLUID
    IF(PRESENT(Troubled_Cell_W)) THEN
      IF (Troubled_Cell_W(ii,jj) .NE. 0) THEN !*PWC
        CYCLE
      END IF
    END IF
#else
PRINT*, "Impossible to use PWCINTROUBLEDCELLS"
PRINT*, "Disable this flag"
PRINT*, "You should not be here"
STOP
#endif
#endif

    SELECT CASE (ABS(ReconstructionUsed))
      CASE(1)
        !*Do nothing. Keep BY_Vol = 0.0 
        !*Because the reconstruction is PWC and the derivatives disappear
        !*NB: Consistency is ensured by the surface contribution of PC
    
      CASE(10,2,20,21,22,23,24,25,26,27,28,29) !*MUSCL, k-MUSCL, MUSCL_CO, MUSCL_VL, MUSCL_M, MUSCL_VA MUSCL_SBM

        VtempMsupport(:,nGPs)=WM(1:nVar,nGPs,ii,jj)
        VtempPsupport(:,nGPs)=WP(1:nVar,nGPs,ii,jj)

#if MOMENTUMINPRIMITIVEVARIABLES
        VtempMsupport(2:3,nGPs)=VtempMsupport(2:3,nGPs)/VtempMsupport(1,nGPs)
        VtempPsupport(2:3,nGPs)=VtempPsupport(2:3,nGPs)/VtempPsupport(1,nGPs)
#endif

        VtempM(1:nVar,nGPs)=VtempMsupport(1:nVar,nGPs)
        VtempP(1:nVar,nGPs)=VtempPsupport(1:nVar,nGPs)


        !*ro 1; vx 2; vy 3; p  4
        dyro = (VtempP(1,nGPs)-VtempM(1,nGPs))

        dyvx = (VtempP(2,nGPs)-VtempM(2,nGPs))
        dyvy = (VtempP(3,nGPs)-VtempM(3,nGPs))

        dyp  = (VtempP(4,nGPs)-VtempM(4,nGPs))
#ifdef MULTIFLUID
#ifdef PRIMITIVEFORMULATIONLEVELSET
        dyPhi= (VtempP(5,nGPs)-VtempM(5,nGPs))
#endif
#endif
        dyrovy = (VtempP(1,nGPs)*VtempP(3,nGPs)-VtempM(1,nGPs)*VtempM(3,nGPs))


        ro = W(1,ii,jj)
        p  = W(4,ii,jj)
#if defined(MOMENTUMINPRIMITIVEVARIABLES)
        !*In this case I am entering with momentum
        vx = W(2,ii,jj)/W(1,ii,jj)
        vy = W(3,ii,jj)/W(1,ii,jj)
#else
        vx = W(2,ii,jj)
        vy = W(3,ii,jj)
#endif

        BVol_in_qp(1,nGPs,nGPs,ii,jj) =  0.



#ifndef MOMENTUMINPRIMITIVEVARIABLES
!*START not MOMENTUMINPRIMITIVEVARIABLES
BVol_in_qp(2,nGPs,nGPs,ii,jj) = -         vy * dyvx   !*NB: vy is already the average
#ifdef IMEX
  !*START IMEX

#ifdef IMEXL2FULLYUPWINDED
    !*START IMEXL2FULLYUPWINDED
    BVol_in_qp(3,nGPs,nGPs,ii,jj) = -     1.0/ro * dyp
#else
    !*ELSE IMEXL2FULLYUPWINDED
    BVol_in_qp(3,nGPs,nGPs,ii,jj) = -     (max_ro-ro)/(ro*max_ro) * dyp
#endif
    !*END IMEXL2FULLYUPWINDED

#else
  !*ELSE IMEX

#if(1==1)
        BVol_in_qp(3,nGPs,nGPs,ii,jj) = -     1.0/ro * dyp                          !*What I suggest to use 
#else
        BVol_in_qp(3,nGPs,nGPs,ii,jj) = -0.5*(1.0/VtempP(1,nGPs)+1.0/VtempM(1,nGPs)) * dyp  !*What Alex initially proposed
#endif

#endif
  !*END IMEX


#ifdef RELAXATION
  !*START RELAXATION
        BVol_in_qp(3,nGPs,nGPs,ii,jj) = BVol_in_qp(3,nGPs,nGPs,ii,jj)/EPS_LM**2 !*In case I divide by EPS_LM**2
#endif
  !*END RELAXATION

#endif
!*END not MOMENTUMINPRIMITIVEVARIABLES


        !*=======================
        !*SPLITTING
        !*=======================
#ifdef MULTIFLUID
roPhi=W(nVar,ii,jj)
Gmm=Get_Gamma(roPhi)
pinf=Get_Pressure_Infinity(roPhi)
#endif

#if defined(ALTERNATIVEFORMULATIONPRESSURE)
!*START ALTERNATIVEFORMULATIONPRESSURE

        ! PRINT*, "This is a consistent formulation."
        ! PRINT*, "It works but it is not suitable for the splitting."
        ! STOP

#ifdef MULTIFLUID
  !*START MULTIFLUID
  BVol_in_qp(4,nGPs,nGPs,ii,jj) = -( (Gmm-1.0)*p + Gmm*pinf ) *dyvy
#else
  !*ELSE MULTIFLUID
  BVol_in_qp(4,nGPs,nGPs,ii,jj) = -( (Gmm-1.0)*p ) *dyvy
#endif
  !*END MULTIFLUID

#ifdef IMEX
  !*START IMEX
  PRINT*, "In particular, for IMEX one really can't do the splitting in this case."
  STOP
#endif
  !*END IMEX

#elif defined(PRESSUREFORMULATIONFORIMEXMOMENTUM)
!*ELSE ALTERNATIVEFORMULATIONPRESSURE

#ifdef IMEXMOMENTUM
  !*START IMEXMOMENTUM
  BVol_in_qp(4,nGPs,nGPs,ii,jj) = Gmm*p/ro*vy * dyro &
                                & -vy *dyp
#else
  !*ELSE IMEXMOMENTUM
  BVol_in_qp(4,nGPs,nGPs,ii,jj) = Gmm*p/ro*vy * dyro &
                                & -Gmm*p/ro* dyrovy &
                                & -vy *dyp
#endif
  !*END IMEXMOMENTUM

#else
!*ELSE ALTERNATIVEFORMULATIONPRESSURE

#ifdef IMEX
  !*START IMEX

#ifdef IMEXL2FULLYUPWINDED
    !*START IMEXL2FULLYUPWINDED
    BVol_in_qp(4,nGPs,nGPs,ii,jj) = - Gmm * p * dyvy &
                                  & -       vy * dyp
#else
    !*ELSE IMEXL2FULLYUPWINDED
    BVol_in_qp(4,nGPs,nGPs,ii,jj) = - Gmm * (p-min_p) * dyvy & !*This must be splitted  
                                  & -       vy * dyp            !*NB: This other one not   
#endif
    !*END IMEXL2FULLYUPWINDED

#else
  !*ELSE IMEX

#ifdef MULTIFLUID
    !*START MULTIFLUID
    BVol_in_qp(4,nGPs,nGPs,ii,jj) = - Gmm * (p+pinf) * dyvy &
                                  & -       vy * dyp
#else
    !*ELSE MULTIFLUID
    BVol_in_qp(4,nGPs,nGPs,ii,jj) = - Gmm * p * dyvy &
                                  & -       vy * dyp
#endif
    !*END MULTIFLUID

#ifdef MULTIFLUID
    !*START MULTIFLUID

#ifndef PRIMITIVEFORMULATIONLEVELSET
      !*START PRIMITIVEFORMULATIONLEVELSET
      BVol_in_qp(nVar,nGPs,nGPs,ii,jj) = 0.0
#else
      !*ELSE PRIMITIVEFORMULATIONLEVELSET
      BVol_in_qp(nVar,nGPs,nGPs,ii,jj) = -vy * dyPhi
#endif
      !*END PRIMITIVEFORMULATIONLEVELSET

#endif
    !*END MULTIFLUID

#endif
  !*END IMEX

#endif
!*END ALTERNATIVEFORMULATIONPRESSURE


        !*--------------------------
        !*SAFETY CHECK
        !*--------------------------
        ! IF (ABS(W(3,ii,jj)-0.5*(VtempP(3,nGPs)+VtempM(3,nGPs))) .GT. 1e-15) THEN
        !   PRINT*, "Problem in PC volume Y contribution"
        !   PRINT*, W(3,ii,jj)-0.5*(VtempP(3,nGPs)+VtempM(3,nGPs))
        !   STOP
        ! END IF
        !*--------------------------

      CASE DEFAULT
        ErrorMessage = "Reconstruction not implemented in BVolY_W_INPUT"
        WRITE(*,*) ErrorMessage
        STOP
    END SELECT
  END DO
END DO



DO jj=lbYs,lbYe
  DO ii=lbXs,lbXe
      DO iGP=1,nGPs
        DO jGP=1,nGPs
          BY_Vol(1:nVar,ii,jj) = BY_Vol(1:nVar,ii,jj) + WeightsGP(iGP,jGP) * BVol_in_qp(1:nVar,iGP,jGP,ii,jj) 
        END DO 
      END DO 
    END DO
  END DO


!-------------------------------------------------------------------------------!
END SUBROUTINE BVolY_W_INPUT
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE UpdateTimeDerivative_W_INPUT(Wt,vbXs,vbXe,vbYs,vbYe)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: FX
USE MOD_FiniteVolume2D_vars,ONLY: FY
USE MOD_FiniteVolume2D_vars,ONLY: BX_Sur
USE MOD_FiniteVolume2D_vars,ONLY: BY_Sur
USE MOD_FiniteVolume2D_vars,ONLY: BX_Vol
USE MOD_FiniteVolume2D_vars,ONLY: BY_Vol
USE MOD_FiniteVolume2D_vars,ONLY: S
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,vbXs:vbXe,vbYs:vbYe), INTENT(INOUT) :: Wt
INTEGER,                                     INTENT(IN)    :: vbXs,vbXe,vbYs,vbYe !*Bounds for the vector
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj
!-------------------------------------------------------------------------------!


DO jj=vbYs,vbYe
  DO ii=vbXs,vbXe
    Wt(1:nVar,ii,jj) =  S(1:nVar,ii,jj) &
                        &  - (FX(1:nVar,ii+0,jj+0)-FX(1:nVar,ii-1,jj+0)-BX_Vol(1:nVar,ii,jj)-BX_Sur(1:nVar,ii-1,jj+0,2)+BX_Sur(1:nVar,ii+0,jj+0,1))/Mesh_DX(1) &
                        &  - (FY(1:nVar,ii+0,jj+0)-FY(1:nVar,ii+0,jj-1)-BY_Vol(1:nVar,ii,jj)-BY_Sur(1:nVar,ii+0,jj-1,2)+BY_Sur(1:nVar,ii+0,jj+0,1))/Mesh_DX(2)
  END DO !ii
END DO !jj


!-------------------------------------------------------------------------------!
END SUBROUTINE UpdateTimeDerivative_W_INPUT
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE AF_PostProcessing_Subroutine(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: AF_PostProcessing
USE MOD_Equation           ,ONLY: BoundaryConditions
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

SELECT CASE (AF_PostProcessing)
  CASE(0)
    !*Do nothing
  CASE(1)
    !*First X, then Y



#ifdef MULTIFLUID
    CALL BoundaryConditions(t)
    CALL Identify_Interface_Cells()
#endif

    CALL BoundaryConditions(t)
    CALL AF_PostProcessing_MINMOD_X()
    CALL BoundaryConditions(t) !*NB: If one does not applies BCs again, conservation is spoiled (clearly).
    CALL AF_PostProcessing_MINMOD_Y()

#ifdef MULTIFLUID
    CALL BoundaryConditions(t)
    CALL Replace_Conserved_At_Interface()
#endif


  CASE(2)
    !*First Y, then X

#ifdef MULTIFLUID
    CALL BoundaryConditions(t)
    CALL Identify_Interface_Cells()
#endif

    CALL BoundaryConditions(t)
    CALL AF_PostProcessing_MINMOD_Y()
    CALL BoundaryConditions(t) !*NB: If one does not applies BCs again, conservation is spoiled (clearly).
    CALL AF_PostProcessing_MINMOD_X()


#ifdef MULTIFLUID
    CALL BoundaryConditions(t)
    CALL Replace_Conserved_At_Interface()
#endif

#ifndef MULTIFLUID
  !*===========================
  !*THE FOLLOWING DEPEND ON EPS
  !*===========================
  CASE(3)
    !*First X, then Y
    PRINT*, "This post-processing is not conservative for relaxation"
    PRINT*, AF_PostProcessing
    STOP
    CALL BoundaryConditions(t)
    CALL AF_PostProcessing_MINMOD_X_DEPENDING_ON_EPS()
    CALL BoundaryConditions(t) !*NB: If one does not applies BCs again, conservation is spoiled (clearly).
    CALL AF_PostProcessing_MINMOD_Y_DEPENDING_ON_EPS()
  CASE(4)
    !*First Y, then X
    PRINT*, "This post-processing is not conservative for relaxation"
    PRINT*, AF_PostProcessing
    STOP
    CALL BoundaryConditions(t)
    CALL AF_PostProcessing_MINMOD_Y_DEPENDING_ON_EPS()
    CALL BoundaryConditions(t) !*NB: If one does not applies BCs again, conservation is spoiled (clearly).
    CALL AF_PostProcessing_MINMOD_X_DEPENDING_ON_EPS()

  CASE(5)
    !*First X, then Y
    CALL BoundaryConditions(t)
    CALL AF_PostProcessing_MINMOD_X_DEPENDING_ON_EPS_CONSERVATIVE()
    CALL BoundaryConditions(t) !*NB: If one does not applies BCs again, conservation is spoiled (clearly).
    CALL AF_PostProcessing_MINMOD_Y_DEPENDING_ON_EPS_CONSERVATIVE()
  CASE(6)
    !*First Y, then X
    CALL BoundaryConditions(t)
    CALL AF_PostProcessing_MINMOD_Y_DEPENDING_ON_EPS_CONSERVATIVE()
    CALL BoundaryConditions(t) !*NB: If one does not applies BCs again, conservation is spoiled (clearly).
    CALL AF_PostProcessing_MINMOD_X_DEPENDING_ON_EPS_CONSERVATIVE()
#endif

  CASE(-1)
    !*First X, then Y



#ifdef MULTIFLUID
    CALL BoundaryConditions(t)
    CALL Identify_Interface_Cells()
#endif

    CALL BoundaryConditions(t)
    CALL AF_PostProcessing_MINMOD_X_Primitive_Based()
    CALL BoundaryConditions(t) !*NB: If one does not applies BCs again, conservation is spoiled (clearly).
    CALL AF_PostProcessing_MINMOD_Y_Primitive_Based()

#ifdef MULTIFLUID
    CALL BoundaryConditions(t)
    CALL Replace_Conserved_At_Interface()
#endif


  CASE(-2)
    !*First Y, then X

#ifdef MULTIFLUID
    CALL BoundaryConditions(t)
    CALL Identify_Interface_Cells()
#endif

    CALL BoundaryConditions(t)
    CALL AF_PostProcessing_MINMOD_Y_Primitive_Based()
    CALL BoundaryConditions(t) !*NB: If one does not applies BCs again, conservation is spoiled (clearly).
    CALL AF_PostProcessing_MINMOD_X_Primitive_Based()


#ifdef MULTIFLUID
    CALL BoundaryConditions(t)
    CALL Replace_Conserved_At_Interface()
#endif

  CASE(11)
    !*First X, then Y



#ifdef MULTIFLUID
    CALL BoundaryConditions(t)
    CALL Identify_Interface_Cells()
#endif

    CALL BoundaryConditions(t)
    CALL AF_PostProcessing_Potentially_HO_X()
    CALL BoundaryConditions(t) !*NB: If one does not applies BCs again, conservation is spoiled (clearly).
    CALL AF_PostProcessing_Potentially_HO_Y()

#ifdef MULTIFLUID
    CALL BoundaryConditions(t)
    CALL Replace_Conserved_At_Interface()
#endif


  CASE(12)
    !*First Y, then X

#ifdef MULTIFLUID
    CALL BoundaryConditions(t)
    CALL Identify_Interface_Cells()
#endif

    CALL BoundaryConditions(t)
    CALL AF_PostProcessing_Potentially_HO_Y()
    CALL BoundaryConditions(t) !*NB: If one does not applies BCs again, conservation is spoiled (clearly).
    CALL AF_PostProcessing_Potentially_HO_X()


#ifdef MULTIFLUID
    CALL BoundaryConditions(t)
    CALL Replace_Conserved_At_Interface()
#endif


  CASE DEFAULT
    ErrorMessage = "AF_PostProcessing not implemented in AF_PostProcessing_Subroutine"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!-------------------------------------------------------------------------------!
END SUBROUTINE AF_PostProcessing_Subroutine
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE AF_PostProcessing_MINMOD_X()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_Equation,           ONLY: PrimToCons
USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Reconstruction,     ONLY: AF_MUSCL
#ifdef MULTIFLUID
USE MOD_FiniteVolume2D_vars,ONLY: Which_Fluid_In_Cell_U
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_U
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_W_X
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_W_Y
USE MOD_Reconstruction,     ONLY: SIDED_AF_MUSCL
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: QtempM(1:nVar)
REAL               :: QtempP(1:nVar)
REAL               :: Qtemp(1:nVar), Vtemp(1:nVar)
REAL               :: Qaverage(1:nVar)
INTEGER            :: final_index
INTEGER            :: fluid_C, fluid_M, fluid_P
CHARACTER          :: which_side
!-------------------------------------------------------------------------------!

#ifdef MULTIFLUID

#ifdef LEVELSETINU
  final_index=nVar
#else
  final_index=nVar-1
#endif

#else

final_index=nVar

#endif

WM=0.0
WP=0.0

SELECT CASE (ABS(Reconstruction))
  CASE(1,10,2,20,21,22,23,24,25,26,27,28,29)

    !*Reconstruction of the limited boundary values

    DO jj=1,nElemsY
      DO ii=0,nElemsX+1 !*NB: Including the ghost cells 0 and nElemsX+1

        !*NB: Limiting taking place in Conserved variables
        CALL PrimToCons(W_X(1:nVar,ii  ,jj),QtempM(1:nVar))
        CALL PrimToCons(W_X(1:nVar,ii+1,jj),QtempP(1:nVar))

        ! PRINT*, ii,jj
        ! PRINT*, W_X(1:nVar,ii  ,jj),W_X(1:nVar,ii+1,jj)
        ! PRINT*, QtempM,             QtempP
        CALL PrimToCons(0.5*(W_X(1:nVar,ii,jj)+W_X(1:nVar,ii+1,jj)),Qaverage(1:nVar))




#ifdef SIDEDPOSTPROCESSING
#ifdef MULTIFLUID
        fluid_C=Which_Fluid_In_Cell_U(ii,jj)
        fluid_M=Which_Fluid_In_Cell_U(ii-1,jj)
        fluid_P=Which_Fluid_In_Cell_U(ii+1,jj)

        IF ((fluid_C .EQ. fluid_M) .AND. (fluid_C .EQ. fluid_P)) THEN
          which_side="C"
        ELSE IF ( (fluid_C .EQ. fluid_M) .AND. (fluid_C .NE. fluid_P) ) THEN
          which_side="M"
        ELSE IF ( (fluid_C .NE. fluid_M) .AND. (fluid_C .EQ. fluid_P) ) THEN
          which_side="P"
        ELSE
          which_side="R"        
        END IF

        CALL SIDED_AF_MUSCL(&
                  U(1:nVar  ,ii-1:ii+1  ,jj), &
                  QtempM(1:nVar),&
                  QtempP(1:nVar),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(1),which_side)
#else
        PRINT*, "Sided PP not existing without multilfuid flag"
        STOP
#endif
#else
        CALL AF_MUSCL(&
                  U(1:nVar  ,ii-1:ii+1  ,jj), &
                  QtempM(1:nVar),&
                  QtempP(1:nVar),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(1),Qaverage(1:nVar))
#endif

        ! PRINT*, ii, jj, QtempM(1:nVar)-WM(1:nVar,1,ii,jj)
        ! PRINT*, ii, jj, QtempP(1:nVar)-WP(1:nVar,1,ii,jj)

      END DO
    END DO

    !*Replacement of the boundary values

    DO jj=1,nElemsY
      DO ii=1,nElemsX+1 !*NB:+1 in X direction

        Qtemp(1:nVar)=0.5*(WP(1:nVar,1,ii-1,jj)+WM(1:nVar,1,ii  ,jj))
        !*+ from the previous cell
        !*- from the actual   cell

        CALL ConsToPrim(Qtemp(1:nVar), Vtemp(1:nVar)) 

#ifdef MULTIFLUID
#ifdef SIDEDPOSTPROCESSING
        IF (Which_Fluid_In_Cell_U(ii-1,jj) .NE. Which_Fluid_In_Cell_U(ii,jj)) THEN
          CYCLE
        END IF
#else
        !*If interface do not touch the primitive variables
        IF (Troubled_Cell_W_X(ii,jj) .NE. 0) THEN
          !*-----------------------
          !*SAFETY CHECK
          !*-----------------------
          ! IF ( (Troubled_Cell_U(ii-1,jj) .EQ. 0) .AND. (Troubled_Cell_U(ii,jj) .EQ. 0) ) THEN
          !   PRINT*, ii, jj
          !   PRINT*, Troubled_Cell_W_X(ii,jj)
          !   PRINT*, Troubled_Cell_U(ii-1,jj), Troubled_Cell_U(ii,jj)
          !   STOP
          ! END IF 
          !*-----------------------

          CYCLE
        END IF
#endif

        !*-----------------------
        !*SAFETY CHECK
        !*-----------------------
        ! IF (Troubled_Cell_W_X(ii,jj) .NE. 1) THEN
        !   IF ( (Troubled_Cell_U(ii-1,jj) .NE. 0) .OR. (Troubled_Cell_U(ii,jj) .NE. 0) ) THEN
        !     PRINT*, ii, jj
        !     PRINT*, Troubled_Cell_W_X(ii,jj)
        !     PRINT*, Troubled_Cell_U(ii-1,jj), Troubled_Cell_U(ii,jj)
        !     STOP
        !   END IF 
        ! END IF 
        !*-----------------------
#endif

        W_X(1:final_index,ii,jj)=Vtemp(1:final_index)

      END DO
    END DO

    !*Replacement of the conserved values
    DO jj=1,nElemsY
      DO ii=1,nElemsX

        CALL PrimToCons(W_X(1:nVar,ii  ,jj),QtempM(1:nVar))
        CALL PrimToCons(W_X(1:nVar,ii+1,jj),QtempP(1:nVar))

#ifdef MULTIFLUID
        !*If interface do not touch the conserved variables
        IF ( (Troubled_Cell_U(ii,jj) .NE. 0) ) THEN
          ! U(2:3,ii,jj)=0.5*(QtempM(2:3)+QtempP(2:3)) !*DO IT ONLY ON VELOCITY
          CYCLE
        END IF
#endif

        U(1:final_index,ii,jj)=0.5*(QtempM(1:final_index)+QtempP(1:final_index))

      END DO
    END DO

  CASE DEFAULT
    ErrorMessage = "AF_PostProcessing_MINMOD_X not implemented for such reconstruction (order higher than 2)"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE AF_PostProcessing_MINMOD_X
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE AF_PostProcessing_MINMOD_Y()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_Equation,           ONLY: PrimToCons
USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Reconstruction,     ONLY: AF_MUSCL
#ifdef MULTIFLUID
USE MOD_FiniteVolume2D_vars,ONLY: Which_Fluid_In_Cell_U
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_U
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_W_X
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_W_Y
USE MOD_Reconstruction,     ONLY: SIDED_AF_MUSCL

#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: QtempM(1:nVar)
REAL               :: QtempP(1:nVar)
REAL               :: Qtemp(1:nVar), Vtemp(1:nVar)
REAL               :: Qaverage(1:nVar)
INTEGER            :: final_index
INTEGER            :: fluid_C, fluid_M, fluid_P
CHARACTER          :: which_side
!-------------------------------------------------------------------------------!

#ifdef MULTIFLUID

#ifdef LEVELSETINU
  final_index=nVar
#else
  final_index=nVar-1
#endif

#else

final_index=nVar

#endif


WM=0.0
WP=0.0

SELECT CASE (ABS(Reconstruction))
  CASE(1,10,2,20,21,22,23,24,25,26,27,28,29)

    !*Reconstruction of the limited boundary values

    DO ii=1,nElemsX
      DO jj=0,nElemsY+1 !*NB: Including the ghost cells 0 and nElemsY+1

        !*NB: Limiting taking place in Conserved variables
        CALL PrimToCons(W_Y(1:nVar,ii,jj  ),QtempM(1:nVar))
        CALL PrimToCons(W_Y(1:nVar,ii,jj+1),QtempP(1:nVar))

        CALL PrimToCons(0.5*(W_Y(1:nVar,ii,jj)+W_Y(1:nVar,ii,jj+1)),Qaverage(1:nVar))

#ifdef SIDEDPOSTPROCESSING
#ifdef MULTIFLUID
        fluid_C=Which_Fluid_In_Cell_U(ii,jj)
        fluid_M=Which_Fluid_In_Cell_U(ii,jj-1)
        fluid_P=Which_Fluid_In_Cell_U(ii,jj+1)

        IF ((fluid_C .EQ. fluid_M) .AND. (fluid_C .EQ. fluid_P)) THEN
          which_side="C"
        ELSE IF ( (fluid_C .EQ. fluid_M) .AND. (fluid_C .NE. fluid_P) ) THEN
          which_side="M"
        ELSE IF ( (fluid_C .NE. fluid_M) .AND. (fluid_C .EQ. fluid_P) ) THEN
          which_side="P"
        ELSE
          which_side="R"        
        END IF

        CALL SIDED_AF_MUSCL(&
                  U(1:nVar  ,ii  ,jj-1:jj+1), &
                  QtempM(1:nVar),&
                  QtempP(1:nVar),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(2),which_side)
#else
        PRINT*, "Sided PP not existing without multilfuid flag"
        STOP
#endif
#else
        CALL AF_MUSCL(&
                  U(1:nVar  ,ii  ,jj-1:jj+1), &
                  QtempM(1:nVar),&
                  QtempP(1:nVar),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(2),Qaverage(1:nVar))
#endif

      END DO
    END DO

    !*Replacement of the boundary values

    DO ii=1,nElemsX
      DO jj=1,nElemsY+1 !*NB:+1 in Y direction

        Qtemp(1:nVar)=0.5*(WP(1:nVar,1,ii,jj-1)+WM(1:nVar,1,ii,jj))
        !*+ from the previous cell
        !*- from the actual   cell

        CALL ConsToPrim(Qtemp(1:nVar), Vtemp(1:nVar)) 

#ifdef MULTIFLUID
#ifdef SIDEDPOSTPROCESSING
        IF (Which_Fluid_In_Cell_U(ii,jj-1) .NE. Which_Fluid_In_Cell_U(ii,jj)) THEN
          CYCLE
        END IF
#else
        !*If interface do not touch the primitive variables
        IF (Troubled_Cell_W_Y(ii,jj) .NE. 0) THEN
          !*-----------------------
          !*SAFETY CHECK
          !*-----------------------
          ! IF ( (Troubled_Cell_U(ii,jj-1) .EQ. 0) .AND. (Troubled_Cell_U(ii,jj) .EQ. 0) ) THEN
          !   PRINT*, ii, jj
          !   PRINT*, Troubled_Cell_W_Y(ii,jj)
          !   PRINT*, Troubled_Cell_U(ii,jj-1), Troubled_Cell_U(ii,jj)
          !   STOP
          ! END IF 
          !*-----------------------

          CYCLE
        END IF
#endif
        !*-----------------------
        !*SAFETY CHECK
        !*-----------------------
        ! IF (Troubled_Cell_W_Y(ii,jj) .NE. 1) THEN
        !   IF ( (Troubled_Cell_U(ii,jj-1) .NE. 0) .OR. (Troubled_Cell_U(ii,jj) .NE. 0) ) THEN
        !     PRINT*, ii, jj
        !     PRINT*, Troubled_Cell_W_Y(ii,jj)
        !     PRINT*, Troubled_Cell_U(ii,jj-1), Troubled_Cell_U(ii,jj)
        !     STOP
        !   END IF 
        ! END IF 
        !*-----------------------
#endif


        W_Y(1:final_index,ii,jj)=Vtemp(1:final_index)

      END DO
    END DO

    !*Replacement of the conserved values
    DO ii=1,nElemsX
      DO jj=1,nElemsY

        CALL PrimToCons(W_Y(1:nVar,ii,jj  ),QtempM(1:nVar))
        CALL PrimToCons(W_Y(1:nVar,ii,jj+1),QtempP(1:nVar))

#ifdef MULTIFLUID
        !*If interface do not touch the conserved variables
        IF ( (Troubled_Cell_U(ii,jj) .NE. 0) ) THEN
          ! U(2:3,ii,jj)=0.5*(QtempM(2:3)+QtempP(2:3)) !*DO IT ONLY ON VELOCITY
          CYCLE
        END IF
#endif

        U(1:final_index,ii,jj)=0.5*(QtempM(1:final_index)+QtempP(1:final_index))


      END DO
    END DO

  CASE DEFAULT
    ErrorMessage = "AF_PostProcessing_MINMOD_Y not implemented for such reconstruction (order higher than 2)"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE AF_PostProcessing_MINMOD_Y
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE AF_PostProcessing_MINMOD_X_DEPENDING_ON_EPS()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_Equation,           ONLY: PrimToCons
USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Reconstruction,     ONLY: AF_MUSCL
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: QtempM(1:nVar)
REAL               :: QtempP(1:nVar)
REAL               :: Qtemp(1:nVar), Vtemp(1:nVar)
!-------------------------------------------------------------------------------!

WM=0.0
WP=0.0

SELECT CASE (ABS(Reconstruction))
  CASE(1,10,2,20,21,22,23,24,25,26,27,28,29)

    !*Reconstruction of the limited boundary values

    DO jj=1,nElemsY
      DO ii=0,nElemsX+1 !*NB: Including the ghost cells 0 and nElemsX+1

        !*NB: Limiting taking place in Conserved variables
        CALL PrimToCons(W_X(1:nVar,ii  ,jj),QtempM(1:nVar))
        CALL PrimToCons(W_X(1:nVar,ii+1,jj),QtempP(1:nVar))

        ! PRINT*, ii,jj
        ! PRINT*, W_X(1:nVar,ii  ,jj),W_X(1:nVar,ii+1,jj)
        ! PRINT*, QtempM,             QtempP

        CALL AF_MUSCL(&
                  U(1:nVar  ,ii-1:ii+1  ,jj), &
                  QtempM(1:nVar),&
                  QtempP(1:nVar),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(1))

        ! PRINT*, ii, jj, QtempM(1:nVar)-WM(1:nVar,1,ii,jj)
        ! PRINT*, ii, jj, QtempP(1:nVar)-WP(1:nVar,1,ii,jj)

      END DO
    END DO

    !*Replacement of the boundary values

    DO jj=1,nElemsY
      DO ii=1,nElemsX+1 !*NB:+1 in X direction

        Qtemp(1:nVar)=0.5*(WP(1:nVar,1,ii-1,jj)+WM(1:nVar,1,ii  ,jj))
        !*+ from the previous cell
        !*- from the actual   cell

        CALL ConsToPrim(Qtemp(1:nVar), Vtemp(1:nVar)) 
        
#ifndef RELAXATION
#ifdef MULTIFLUID
        W_X(1:nVar-1,ii,jj)=Vtemp(1:nVar-1)
#else
        W_X(1:nVar,ii,jj)=Vtemp(1:nVar)
#endif
#else
        W_X(1:nVar,ii,jj)=Vtemp(1:nVar)*EPS_LM**2+(1.0-EPS_LM**2)*W_X(1:nVar,ii,jj)
#endif

      END DO
    END DO

    !*Replacement of the conserved values
    DO jj=1,nElemsY
      DO ii=1,nElemsX

        CALL PrimToCons(W_X(1:nVar,ii  ,jj),QtempM(1:nVar))
        CALL PrimToCons(W_X(1:nVar,ii+1,jj),QtempP(1:nVar))

#ifndef RELAXATION
#ifdef MULTIFLUID
        U(1:nVar-1,ii,jj)=0.5*(QtempM(1:nVar-1)+QtempP(1:nVar-1))
#else
        U(1:nVar,ii,jj)=0.5*(QtempM(1:nVar)+QtempP(1:nVar))
#endif

#else
        U(1:nVar,ii,jj)=0.5*(QtempM(1:nVar)+QtempP(1:nVar))*EPS_LM**2+(1.0-EPS_LM**2)*U(1:nVar,ii,jj)
#endif
      END DO
    END DO

  CASE DEFAULT
    ErrorMessage = "AF_PostProcessing_MINMOD_X_DEPENDING_ON_EPS not implemented for such reconstruction (order higher than 2)"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE AF_PostProcessing_MINMOD_X_DEPENDING_ON_EPS
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE AF_PostProcessing_MINMOD_Y_DEPENDING_ON_EPS()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_Equation,           ONLY: PrimToCons
USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Reconstruction,     ONLY: AF_MUSCL
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: QtempM(1:nVar)
REAL               :: QtempP(1:nVar)
REAL               :: Qtemp(1:nVar), Vtemp(1:nVar)
!-------------------------------------------------------------------------------!

WM=0.0
WP=0.0

SELECT CASE (ABS(Reconstruction))
  CASE(1,10,2,20,21,22,23,24,25,26,27,28,29)

    !*Reconstruction of the limited boundary values

    DO ii=1,nElemsX
      DO jj=0,nElemsY+1 !*NB: Including the ghost cells 0 and nElemsY+1

        !*NB: Limiting taking place in Conserved variables
        CALL PrimToCons(W_Y(1:nVar,ii,jj  ),QtempM(1:nVar))
        CALL PrimToCons(W_Y(1:nVar,ii,jj+1),QtempP(1:nVar))


        CALL AF_MUSCL(&
                  U(1:nVar  ,ii  ,jj-1:jj+1), &
                  QtempM(1:nVar),&
                  QtempP(1:nVar),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(2))


      END DO
    END DO

    !*Replacement of the boundary values

    DO ii=1,nElemsX
      DO jj=1,nElemsY+1 !*NB:+1 in Y direction

        Qtemp(1:nVar)=0.5*(WP(1:nVar,1,ii,jj-1)+WM(1:nVar,1,ii,jj))
        !*+ from the previous cell
        !*- from the actual   cell

        CALL ConsToPrim(Qtemp(1:nVar), Vtemp(1:nVar)) 

#ifndef RELAXATION        
#ifdef MULTIFLUID
        W_Y(1:nVar-1,ii,jj)=Vtemp(1:nVar-1)
#else
        W_Y(1:nVar,ii,jj)=Vtemp(1:nVar)
#endif
#else
        W_Y(1:nVar,ii,jj)=Vtemp(1:nVar)*EPS_LM**2+(1.0-EPS_LM**2)*W_Y(1:nVar,ii,jj)
#endif

      END DO
    END DO

    !*Replacement of the conserved values
    DO ii=1,nElemsX
      DO jj=1,nElemsY

        CALL PrimToCons(W_Y(1:nVar,ii,jj  ),QtempM(1:nVar))
        CALL PrimToCons(W_Y(1:nVar,ii,jj+1),QtempP(1:nVar))

#ifndef RELAXATION
#ifdef MULTIFLUID
        U(1:nVar-1,ii,jj)=0.5*(QtempM(1:nVar-1)+QtempP(1:nVar-1))
#else
        U(1:nVar,ii,jj)=0.5*(QtempM(1:nVar)+QtempP(1:nVar))
#endif
#else
        U(1:nVar,ii,jj)=0.5*(QtempM(1:nVar)+QtempP(1:nVar))*EPS_LM**2+(1.0-EPS_LM**2)*U(1:nVar,ii,jj)
#endif

      END DO
    END DO

  CASE DEFAULT
    ErrorMessage = "AF_PostProcessing_MINMOD_Y_DEPENDING_ON_EPS not implemented for such reconstruction (order higher than 2)"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE AF_PostProcessing_MINMOD_Y_DEPENDING_ON_EPS
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE AF_PostProcessing_MINMOD_X_DEPENDING_ON_EPS_CONSERVATIVE()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_Equation,           ONLY: PrimToCons
USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Reconstruction,     ONLY: AF_MUSCL
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: QtempM(1:nVar)
REAL               :: QtempP(1:nVar)
REAL               :: Qtemp(1:nVar), Vtemp(1:nVar)
REAL               :: Q_auxiliary_for_conservation(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
!-------------------------------------------------------------------------------!

WM=0.0
WP=0.0

SELECT CASE (ABS(Reconstruction))
  CASE(1,10,2,20,21,22,23,24,25,26,27,28,29)

    !*Reconstruction of the limited boundary values

    DO jj=1,nElemsY
      DO ii=0,nElemsX+1 !*NB: Including the ghost cells 0 and nElemsX+1

        !*NB: Limiting taking place in Conserved variables
        CALL PrimToCons(W_X(1:nVar,ii  ,jj),QtempM(1:nVar))
        CALL PrimToCons(W_X(1:nVar,ii+1,jj),QtempP(1:nVar))

        ! PRINT*, ii,jj
        ! PRINT*, W_X(1:nVar,ii  ,jj),W_X(1:nVar,ii+1,jj)
        ! PRINT*, QtempM,             QtempP

        CALL AF_MUSCL(&
                  U(1:nVar  ,ii-1:ii+1  ,jj), &
                  QtempM(1:nVar),&
                  QtempP(1:nVar),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(1))

        ! PRINT*, ii, jj, QtempM(1:nVar)-WM(1:nVar,1,ii,jj)
        ! PRINT*, ii, jj, QtempP(1:nVar)-WP(1:nVar,1,ii,jj)

      END DO
    END DO

    !*Replacement of the boundary values

    DO jj=1,nElemsY
      DO ii=1,nElemsX+1 !*NB:+1 in X direction

        Qtemp(1:nVar)=0.5*(WP(1:nVar,1,ii-1,jj)+WM(1:nVar,1,ii  ,jj))
        !*+ from the previous cell
        !*- from the actual   cell

        CALL ConsToPrim(Qtemp(1:nVar), Vtemp(1:nVar)) 

        !*============================
        !*W_X really needs to be postprocessed with this blending to get AP
        !*============================
#ifndef RELAXATION
#ifdef MULTIFLUID
        W_X(1:nVar-1,ii,jj)=Vtemp(1:nVar-1)
#else
        W_X(1:nVar,ii,jj)=Vtemp(1:nVar)
#endif
#else
        W_X(1:nVar,ii,jj)=Vtemp(1:nVar)*EPS_LM**2+(1.0-EPS_LM**2)*W_X(1:nVar,ii,jj)
#endif

        !*============================
        !*However, we also need the unblended W_X for ensuring conservation
        !*============================
        Q_auxiliary_for_conservation(1:nVar,ii,jj)=Qtemp(1:nVar)

      END DO
    END DO

    !*Replacement of the conserved values
    DO jj=1,nElemsY
      DO ii=1,nElemsX

        QtempM(1:nVar) = Q_auxiliary_for_conservation(1:nVar,ii  ,jj)
        QtempP(1:nVar) = Q_auxiliary_for_conservation(1:nVar,ii+1,jj)

#ifndef RELAXATION
#ifdef MULTIFLUID
        U(1:nVar-1,ii,jj)=0.5*(QtempM(1:nVar-1)+QtempP(1:nVar-1))
#else
        U(1:nVar,ii,jj)=0.5*(QtempM(1:nVar)+QtempP(1:nVar))
#endif
#else
        U(1:nVar,ii,jj)=0.5*(QtempM(1:nVar)+QtempP(1:nVar))*EPS_LM**2+(1.0-EPS_LM**2)*U(1:nVar,ii,jj)
#endif
      END DO
    END DO

  CASE DEFAULT
    ErrorMessage = "AF_PostProcessing_MINMOD_X_DEPENDING_ON_EPS_CONSERVATIVE not implemented for such reconstruction (order higher than 2)"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE AF_PostProcessing_MINMOD_X_DEPENDING_ON_EPS_CONSERVATIVE
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE AF_PostProcessing_MINMOD_Y_DEPENDING_ON_EPS_CONSERVATIVE()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_Equation,           ONLY: PrimToCons
USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Reconstruction,     ONLY: AF_MUSCL
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: QtempM(1:nVar)
REAL               :: QtempP(1:nVar)
REAL               :: Qtemp(1:nVar), Vtemp(1:nVar)
REAL               :: Q_auxiliary_for_conservation(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
!-------------------------------------------------------------------------------!

WM=0.0
WP=0.0

SELECT CASE (ABS(Reconstruction))
  CASE(1,10,2,20,21,22,23,24,25,26,27,28,29)

    !*Reconstruction of the limited boundary values

    DO ii=1,nElemsX
      DO jj=0,nElemsY+1 !*NB: Including the ghost cells 0 and nElemsY+1

        !*NB: Limiting taking place in Conserved variables
        CALL PrimToCons(W_Y(1:nVar,ii,jj  ),QtempM(1:nVar))
        CALL PrimToCons(W_Y(1:nVar,ii,jj+1),QtempP(1:nVar))


        CALL AF_MUSCL(&
                  U(1:nVar  ,ii  ,jj-1:jj+1), &
                  QtempM(1:nVar),&
                  QtempP(1:nVar),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(2))


      END DO
    END DO

    !*Replacement of the boundary values

    DO ii=1,nElemsX
      DO jj=1,nElemsY+1 !*NB:+1 in Y direction

        Qtemp(1:nVar)=0.5*(WP(1:nVar,1,ii,jj-1)+WM(1:nVar,1,ii,jj))
        !*+ from the previous cell
        !*- from the actual   cell

        CALL ConsToPrim(Qtemp(1:nVar), Vtemp(1:nVar)) 

        !*============================
        !*W_Y really needs to be postprocessed with this blending to get AP
        !*============================
#ifndef RELAXATION        
#ifdef MULTIFLUID
        W_Y(1:nVar-1,ii,jj)=Vtemp(1:nVar-1)
#else
        W_Y(1:nVar,ii,jj)=Vtemp(1:nVar)
#endif
#else
        W_Y(1:nVar,ii,jj)=Vtemp(1:nVar)*EPS_LM**2+(1.0-EPS_LM**2)*W_Y(1:nVar,ii,jj)
#endif

        !*============================
        !*However, we also need the unblended W_Y for ensuring conservation
        !*============================
        Q_auxiliary_for_conservation(1:nVar,ii,jj)=Qtemp(1:nVar)

      END DO
    END DO

    !*Replacement of the conserved values
    DO ii=1,nElemsX
      DO jj=1,nElemsY

        QtempM(1:nVar) = Q_auxiliary_for_conservation(1:nVar,ii,jj  )
        QtempP(1:nVar) = Q_auxiliary_for_conservation(1:nVar,ii,jj+1)

#ifndef RELAXATION
#ifdef MULTIFLUID
        U(1:nVar-1,ii,jj)=0.5*(QtempM(1:nVar-1)+QtempP(1:nVar-1))
#else
        U(1:nVar,ii,jj)=0.5*(QtempM(1:nVar)+QtempP(1:nVar))
#endif
#else
        U(1:nVar,ii,jj)=0.5*(QtempM(1:nVar)+QtempP(1:nVar))*EPS_LM**2+(1.0-EPS_LM**2)*U(1:nVar,ii,jj)
#endif

      END DO
    END DO

  CASE DEFAULT
    ErrorMessage = "AF_PostProcessing_MINMOD_Y_DEPENDING_ON_EPS_CONSERVATIVE not implemented for such reconstruction (order higher than 2)"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE AF_PostProcessing_MINMOD_Y_DEPENDING_ON_EPS_CONSERVATIVE
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE AF_PostProcessing_MINMOD_X_Primitive_Based()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_Equation,           ONLY: PrimToCons
USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Reconstruction,     ONLY: AF_MUSCL_Primitive_Based
#ifdef MULTIFLUID
USE MOD_FiniteVolume2D_vars,ONLY: Which_Fluid_In_Cell_U
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_U
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_W_X
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_W_Y
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: QtempM(1:nVar)
REAL               :: QtempP(1:nVar)
REAL               :: Qtemp(1:nVar), Vtemp(1:nVar), Uhat(1:nVar)
REAL               :: Qaverage(1:nVar)
REAL               :: Fluctuation(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) !*Same size as U
INTEGER            :: final_index
INTEGER            :: fluid_C, fluid_M, fluid_P
CHARACTER          :: which_side
!-------------------------------------------------------------------------------!

#ifdef MULTIFLUID

#ifdef LEVELSETINU
  final_index=nVar
#else
  final_index=nVar-1
#endif

#else

final_index=nVar

#endif

WM=0.0
WP=0.0
Qaverage=0.0

SELECT CASE (ABS(Reconstruction))
  CASE(1,10,2,20,21,22,23,24,25,26,27,28,29)

    !*Computation of the fluctuations between averages of conserved values and same quantity obtained from primtive variables

    DO jj=1,nElemsY
      DO ii=-1,nElemsX+2 !*NB: Including one layer more with respect to the ghost cells 0 and nElemsX+1

        CALL PrimToCons(W_X(1:nVar,ii  ,jj),QtempM(1:nVar))
        CALL PrimToCons(W_X(1:nVar,ii+1,jj),QtempP(1:nVar))

        ! PRINT*, ii,jj
        ! PRINT*, W_X(1:nVar,ii  ,jj),W_X(1:nVar,ii+1,jj)
        ! PRINT*, QtempM,             QtempP
        Qaverage(1:nVar)=0.5*(QtempM(1:nVar)+QtempP(1:nVar))

        Fluctuation(1:nVar,ii,jj)=U(1:nVar,ii,jj)-Qaverage(1:nVar)

      END DO
    END DO

    !*Linear reconstruction of fluctuation

    DO jj=1,nElemsY
      DO ii=0,nElemsX+1 !*NB: including the ghost cells 0 and nElemsX+1

#ifdef SIDEDPOSTPROCESSING
#ifdef MULTIFLUID
        fluid_C=Which_Fluid_In_Cell_U(ii,jj)
        fluid_M=Which_Fluid_In_Cell_U(ii-1,jj)
        fluid_P=Which_Fluid_In_Cell_U(ii+1,jj)

        IF ((fluid_C .EQ. fluid_M) .AND. (fluid_C .EQ. fluid_P)) THEN
          which_side="C"
        ELSE IF ( (fluid_C .EQ. fluid_M) .AND. (fluid_C .NE. fluid_P) ) THEN
          which_side="M"
        ELSE IF ( (fluid_C .NE. fluid_M) .AND. (fluid_C .EQ. fluid_P) ) THEN
          which_side="P"
        ELSE
          which_side="R"        
        END IF

        ! CALL SIDED_AF_MUSCL_Primitive_Based(&
        !           Fluctuation(1:nVar  ,ii-1:ii+1  ,jj), &
        !           QtempM(1:nVar),&
        !           QtempP(1:nVar),&
        !           WM(1:nVar,1:nGPs,ii,jj),&
        !           WP(1:nVar,1:nGPs,ii,jj),&
        !           MESH_DX(1),which_side)

        PRINT*, "Sided PP possible in theory but not coded yet"

#else
        PRINT*, "Sided PP not existing without multilfuid flag"
        STOP
#endif
#else
        CALL AF_MUSCL_Primitive_Based(&
                  Fluctuation(1:nVar  ,ii-1:ii+1  ,jj), &
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(1))
#endif


      END DO
    END DO

    !*Replacement of the boundary values

    DO jj=1,nElemsY
      DO ii=1,nElemsX+1 !*NB:+1 in X direction

        Vtemp(1:nVar)=0.5*(WP(1:nVar,1,ii-1,jj)+WM(1:nVar,1,ii  ,jj))
        !*+ from the previous cell
        !*- from the actual   cell

        CALL PrimToCons(W_X(1:nVar,ii,jj), Qtemp(1:nVar)) 

        Uhat(1:nVar)=Vtemp(1:nVar)+Qtemp(1:nVar)

#ifdef MULTIFLUID
#ifdef SIDEDPOSTPROCESSING
        IF (Which_Fluid_In_Cell_U(ii-1,jj) .NE. Which_Fluid_In_Cell_U(ii,jj)) THEN
          CYCLE
        END IF
#else
        !*If interface do not touch the primitive variables
        IF (Troubled_Cell_W_X(ii,jj) .NE. 0) THEN
          !*-----------------------
          !*SAFETY CHECK
          !*-----------------------
          ! IF ( (Troubled_Cell_U(ii-1,jj) .EQ. 0) .AND. (Troubled_Cell_U(ii,jj) .EQ. 0) ) THEN
          !   PRINT*, ii, jj
          !   PRINT*, Troubled_Cell_W_X(ii,jj)
          !   PRINT*, Troubled_Cell_U(ii-1,jj), Troubled_Cell_U(ii,jj)
          !   STOP
          ! END IF 
          !*-----------------------

          CYCLE
        END IF
#endif

        !*-----------------------
        !*SAFETY CHECK
        !*-----------------------
        ! IF (Troubled_Cell_W_X(ii,jj) .NE. 1) THEN
        !   IF ( (Troubled_Cell_U(ii-1,jj) .NE. 0) .OR. (Troubled_Cell_U(ii,jj) .NE. 0) ) THEN
        !     PRINT*, ii, jj
        !     PRINT*, Troubled_Cell_W_X(ii,jj)
        !     PRINT*, Troubled_Cell_U(ii-1,jj), Troubled_Cell_U(ii,jj)
        !     STOP
        !   END IF 
        ! END IF 
        !*-----------------------
#endif

        CALL ConsToPrim(Uhat(1:nVar),W_X(1:nVar,ii,jj))

      END DO
    END DO

    !*Replacement of the conserved values
    DO jj=1,nElemsY
      DO ii=1,nElemsX

        CALL PrimToCons(W_X(1:nVar,ii  ,jj),QtempM(1:nVar))
        CALL PrimToCons(W_X(1:nVar,ii+1,jj),QtempP(1:nVar))

#ifdef MULTIFLUID
        !*If interface do not touch the conserved variables
        IF ( (Troubled_Cell_U(ii,jj) .NE. 0) ) THEN
          ! U(2:3,ii,jj)=0.5*(QtempM(2:3)+QtempP(2:3)) !*DO IT ONLY ON VELOCITY
          CYCLE
        END IF
#endif

        U(1:final_index,ii,jj)=0.5*(QtempM(1:final_index)+QtempP(1:final_index))

      END DO
    END DO

  CASE DEFAULT
    ErrorMessage = "AF_PostProcessing_MINMOD_X_Primitive_Based not implemented for such reconstruction (order higher than 2)"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE AF_PostProcessing_MINMOD_X_Primitive_Based
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE AF_PostProcessing_MINMOD_Y_Primitive_Based()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_Equation,           ONLY: PrimToCons
USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Reconstruction,     ONLY: AF_MUSCL_Primitive_Based
#ifdef MULTIFLUID
USE MOD_FiniteVolume2D_vars,ONLY: Which_Fluid_In_Cell_U
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_U
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_W_X
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_W_Y
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: QtempM(1:nVar)
REAL               :: QtempP(1:nVar)
REAL               :: Qtemp(1:nVar), Vtemp(1:nVar), Uhat(1:nVar)
REAL               :: Qaverage(1:nVar)
REAL               :: Fluctuation(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) !*Same size as U
INTEGER            :: final_index
INTEGER            :: fluid_C, fluid_M, fluid_P
CHARACTER          :: which_side
!-------------------------------------------------------------------------------!

#ifdef MULTIFLUID

#ifdef LEVELSETINU
  final_index=nVar
#else
  final_index=nVar-1
#endif

#else

final_index=nVar

#endif


WM=0.0
WP=0.0

SELECT CASE (ABS(Reconstruction))
  CASE(1,10,2,20,21,22,23,24,25,26,27,28,29)

    !*Computation of the fluctuations between averages of conserved values and same quantity obtained from primtive variables

    DO ii=1,nElemsX
      DO jj=-1,nElemsY+2 !*NB: Including one layer more with respect to the ghost cells 0 and nElemsY+1

        CALL PrimToCons(W_Y(1:nVar,ii,jj  ),QtempM(1:nVar))
        CALL PrimToCons(W_Y(1:nVar,ii,jj+1),QtempP(1:nVar))

        Qaverage(1:nVar)=0.5*(QtempM(1:nVar)+QtempP(1:nVar))

        Fluctuation(1:nVar,ii,jj)=U(1:nVar,ii,jj)-Qaverage(1:nVar)

      END DO
    END DO

    !*Linear reconstruction of fluctuation

    DO ii=1,nElemsX
      DO jj=0,nElemsY+1 !*NB: Including the ghost cells 0 and nElemsY+1

#ifdef SIDEDPOSTPROCESSING
#ifdef MULTIFLUID
        fluid_C=Which_Fluid_In_Cell_U(ii,jj)
        fluid_M=Which_Fluid_In_Cell_U(ii,jj-1)
        fluid_P=Which_Fluid_In_Cell_U(ii,jj+1)

        IF ((fluid_C .EQ. fluid_M) .AND. (fluid_C .EQ. fluid_P)) THEN
          which_side="C"
        ELSE IF ( (fluid_C .EQ. fluid_M) .AND. (fluid_C .NE. fluid_P) ) THEN
          which_side="M"
        ELSE IF ( (fluid_C .NE. fluid_M) .AND. (fluid_C .EQ. fluid_P) ) THEN
          which_side="P"
        ELSE
          which_side="R"        
        END IF

        ! CALL SIDED_AF_MUSCL_Primitive_Based(&
        !           U(1:nVar  ,ii  ,jj-1:jj+1), &
        !           QtempM(1:nVar),&
        !           QtempP(1:nVar),&
        !           WM(1:nVar,1:nGPs,ii,jj),&
        !           WP(1:nVar,1:nGPs,ii,jj),&
        !           MESH_DX(2),which_side)

        PRINT*, "Sided PP possible in theory but not coded yet"

#else
        PRINT*, "Sided PP not existing without multilfuid flag"
        STOP
#endif
#else
        CALL AF_MUSCL_Primitive_Based(&
                  Fluctuation(1:nVar  ,ii  ,jj-1:jj+1), &
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(2))
#endif

      END DO
    END DO

    !*Replacement of the boundary values

    DO ii=1,nElemsX
      DO jj=1,nElemsY+1 !*NB:+1 in Y direction



        Vtemp(1:nVar)=0.5*(WP(1:nVar,1,ii,jj-1)+WM(1:nVar,1,ii,jj))
        !*+ from the previous cell
        !*- from the actual   cell

        CALL PrimToCons(W_Y(1:nVar,ii,jj), Qtemp(1:nVar)) 

        Uhat(1:nVar)=Vtemp(1:nVar)+Qtemp(1:nVar)

#ifdef MULTIFLUID
#ifdef SIDEDPOSTPROCESSING
        IF (Which_Fluid_In_Cell_U(ii,jj-1) .NE. Which_Fluid_In_Cell_U(ii,jj)) THEN
          CYCLE
        END IF
#else
        !*If interface do not touch the primitive variables
        IF (Troubled_Cell_W_Y(ii,jj) .NE. 0) THEN
          !*-----------------------
          !*SAFETY CHECK
          !*-----------------------
          ! IF ( (Troubled_Cell_U(ii,jj-1) .EQ. 0) .AND. (Troubled_Cell_U(ii,jj) .EQ. 0) ) THEN
          !   PRINT*, ii, jj
          !   PRINT*, Troubled_Cell_W_Y(ii,jj)
          !   PRINT*, Troubled_Cell_U(ii,jj-1), Troubled_Cell_U(ii,jj)
          !   STOP
          ! END IF 
          !*-----------------------

          CYCLE
        END IF
#endif
        !*-----------------------
        !*SAFETY CHECK
        !*-----------------------
        ! IF (Troubled_Cell_W_Y(ii,jj) .NE. 1) THEN
        !   IF ( (Troubled_Cell_U(ii,jj-1) .NE. 0) .OR. (Troubled_Cell_U(ii,jj) .NE. 0) ) THEN
        !     PRINT*, ii, jj
        !     PRINT*, Troubled_Cell_W_Y(ii,jj)
        !     PRINT*, Troubled_Cell_U(ii,jj-1), Troubled_Cell_U(ii,jj)
        !     STOP
        !   END IF 
        ! END IF 
        !*-----------------------
#endif

        CALL ConsToPrim(Uhat(1:nVar),W_Y(1:nVar,ii,jj))

      END DO
    END DO

    !*Replacement of the conserved values
    DO ii=1,nElemsX
      DO jj=1,nElemsY

        CALL PrimToCons(W_Y(1:nVar,ii,jj  ),QtempM(1:nVar))
        CALL PrimToCons(W_Y(1:nVar,ii,jj+1),QtempP(1:nVar))

#ifdef MULTIFLUID
        !*If interface do not touch the conserved variables
        IF ( (Troubled_Cell_U(ii,jj) .NE. 0) ) THEN
          ! U(2:3,ii,jj)=0.5*(QtempM(2:3)+QtempP(2:3)) !*DO IT ONLY ON VELOCITY
          CYCLE
        END IF
#endif

        U(1:final_index,ii,jj)=0.5*(QtempM(1:final_index)+QtempP(1:final_index))


      END DO
    END DO

  CASE DEFAULT
    ErrorMessage = "AF_PostProcessing_MINMOD_Y_Primitive_Based not implemented for such reconstruction (order higher than 2)"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE AF_PostProcessing_MINMOD_Y_Primitive_Based
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE AF_PostProcessing_Potentially_HO_X()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_Equation,           ONLY: PrimToCons
USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Reconstruction,     ONLY: AF_MUSCL
USE MOD_Reconstruction,     ONLY: MINMOD_3_INPUTS
USE MOD_Reconstruction,     ONLY: SBMMUSCL
#ifdef MULTIFLUID
USE MOD_FiniteVolume2D_vars,ONLY: Which_Fluid_In_Cell_U
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_U
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_W_X
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_W_Y
USE MOD_Reconstruction,     ONLY: SIDED_AF_MUSCL
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj, iVar
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: QtempM(1:nVar)
REAL               :: QtempP(1:nVar)
REAL               :: Qtemp(1:nVar), Vtemp(1:nVar)
REAL               :: Qaverage(1:nVar)
INTEGER            :: final_index
INTEGER            :: fluid_C, fluid_M, fluid_P
CHARACTER          :: which_side
REAL               :: slopes_U(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
REAL               :: slopes_W_X( 1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
REAL               :: Cons_W_X( 1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
REAL               :: Q_vector_dummy( 1:nVar,-1:+1)
REAL               :: theta=1.3
REAL               :: tau=0.5
REAL               :: WM_dummy(1:nVar,1:nGPs)
REAL               :: WP_dummy(1:nVar,1:nGPs)
!-------------------------------------------------------------------------------!

#ifdef MULTIFLUID

#ifdef LEVELSETINU
  final_index=nVar
#else
  final_index=nVar-1
#endif

#else

final_index=nVar

#endif

WM=0.0
WP=0.0
slopes_U  =0.0
slopes_W_X=0.0
Cons_W_X  =0.0

SELECT CASE (ABS(Reconstruction))
  CASE(1,10,2,20,21,22,23,24,25,26,27,28,29)

    !*Reconstruction of the limited boundary values

    DO jj=1,nElemsY
      DO ii=-1,nElemsX+2 

        !*NB: Limiting taking place in Conserved variables
        CALL PrimToCons(W_X(1:nVar,ii  ,jj),QtempM(1:nVar))
        CALL PrimToCons(W_X(1:nVar,ii+1,jj),QtempP(1:nVar))


#ifdef SIDEDPOSTPROCESSING
#ifdef MULTIFLUID
        fluid_C=Which_Fluid_In_Cell_U(ii,jj)
        fluid_M=Which_Fluid_In_Cell_U(ii-1,jj)
        fluid_P=Which_Fluid_In_Cell_U(ii+1,jj)

        IF ((fluid_C .EQ. fluid_M) .AND. (fluid_C .EQ. fluid_P)) THEN
          which_side="C"
        ELSE IF ( (fluid_C .EQ. fluid_M) .AND. (fluid_C .NE. fluid_P) ) THEN
          which_side="M"
        ELSE IF ( (fluid_C .NE. fluid_M) .AND. (fluid_C .EQ. fluid_P) ) THEN
          which_side="P"
        ELSE
          which_side="R"        
        END IF

        ! CALL SIDED_AF_MUSCL(&
        !           U(1:nVar  ,ii-1:ii+1  ,jj), &
        !           QtempM(1:nVar),&
        !           QtempP(1:nVar),&
        !           WM(1:nVar,1:nGPs,ii,jj),&
        !           WP(1:nVar,1:nGPs,ii,jj),&
        !           MESH_DX(1),which_side)

        PRINT*, "Sided PP not coded yet for this PP"
        STOP        
#else
        PRINT*, "Sided PP not existing without multilfuid flag"
        STOP
#endif
#else
        Q_vector_dummy(1:nVar,-1)=QtempM(1:nVar)
        Q_vector_dummy(1:nVar, 0)=U(1:nVar,ii,jj)
        Q_vector_dummy(1:nVar,+1)=QtempP(1:nVar)

#ifdef OVERCOMPRESSIVEPOSTPROCESSING
        CALL SBMMUSCL(Q_vector_dummy(1:nVar,-1:+1),WM_dummy(1:nVar,1:nGPs),WP_dummy(1:nVar,1:nGPs),MESH_DX(1)/2.0,slopes_U(1:nVar,ii,jj),tau_input=tau,theta_input=theta)
#else
        DO iVar=1,nVar
          slopes_U(iVar,ii,jj)=MINMOD_3_INPUTS((Q_vector_dummy(iVar, 0)-Q_vector_dummy(iVar,-1))/(0.5*MESH_DX(1))*theta, &
                                            &  (Q_vector_dummy(iVar,+1)-Q_vector_dummy(iVar,-1))/(    MESH_DX(1)), &
                                            &  (Q_vector_dummy(iVar,+1)-Q_vector_dummy(iVar, 0))/(0.5*MESH_DX(1))*theta &
                                            &  )
        END DO
#endif




#endif


      END DO
    END DO

    !*Replacement of the boundary values

    DO jj=1,nElemsY
      DO ii=0,nElemsX+2 

        Cons_W_X(1:nVar,ii,jj)=0.5*(U(1:nVar,ii-1,jj)+U(1:nVar,ii  ,jj))+MESH_DX(1)/8.0*(slopes_U(1:nVar,ii-1,jj)-slopes_U(1:nVar,ii  ,jj))
        !*from the previous cell
        !*from the actual   cell

#ifdef MULTIFLUID
#ifdef SIDEDPOSTPROCESSING
        IF (Which_Fluid_In_Cell_U(ii-1,jj) .NE. Which_Fluid_In_Cell_U(ii,jj)) THEN
          CYCLE
        END IF
#else
        !*If interface do not touch the primitive variables
        IF (Troubled_Cell_W_X(ii,jj) .NE. 0) THEN
          !*-----------------------
          !*SAFETY CHECK
          !*-----------------------
          ! IF ( (Troubled_Cell_U(ii-1,jj) .EQ. 0) .AND. (Troubled_Cell_U(ii,jj) .EQ. 0) ) THEN
          !   PRINT*, ii, jj
          !   PRINT*, Troubled_Cell_W_X(ii,jj)
          !   PRINT*, Troubled_Cell_U(ii-1,jj), Troubled_Cell_U(ii,jj)
          !   STOP
          ! END IF 
          !*-----------------------

          CYCLE
        END IF
#endif

        !*-----------------------
        !*SAFETY CHECK
        !*-----------------------
        ! IF (Troubled_Cell_W_X(ii,jj) .NE. 1) THEN
        !   IF ( (Troubled_Cell_U(ii-1,jj) .NE. 0) .OR. (Troubled_Cell_U(ii,jj) .NE. 0) ) THEN
        !     PRINT*, ii, jj
        !     PRINT*, Troubled_Cell_W_X(ii,jj)
        !     PRINT*, Troubled_Cell_U(ii-1,jj), Troubled_Cell_U(ii,jj)
        !     STOP
        !   END IF 
        ! END IF 
        !*-----------------------
#endif

        CALL ConsToPrim(Cons_W_X(1:nVar,ii,jj), W_X(1:nVar,ii,jj)) 

      END DO
    END DO

    !*Slopes of the staggered mesh

    DO jj=1,nElemsY
      DO ii=1,nElemsX+1

        Q_vector_dummy(1:nVar,-1)=Cons_W_X(1:nVar,ii-1,jj)
        Q_vector_dummy(1:nVar, 0)=Cons_W_X(1:nVar,ii,jj)
        Q_vector_dummy(1:nVar,+1)=Cons_W_X(1:nVar,ii+1,jj)

#ifdef OVERCOMPRESSIVEPOSTPROCESSING
        CALL SBMMUSCL(Q_vector_dummy(1:nVar,-1:+1),WM_dummy(1:nVar,1:nGPs),WP_dummy(1:nVar,1:nGPs),MESH_DX(1),slopes_W_X(1:nVar,ii,jj),tau_input=tau,theta_input=theta)
#else
        DO iVar=1,nVar
          slopes_W_X(iVar,ii,jj)=MINMOD_3_INPUTS((Q_vector_dummy(iVar, 0)-Q_vector_dummy(iVar,-1))/(    MESH_DX(1))*theta, &
                                              &  (Q_vector_dummy(iVar,+1)-Q_vector_dummy(iVar,-1))/(2.0*MESH_DX(1)), &
                                              &  (Q_vector_dummy(iVar,+1)-Q_vector_dummy(iVar, 0))/(    MESH_DX(1))*theta &
                                              &  )

        END DO
#endif

      END DO
    END DO

    !*Replacement of the conserved values
    DO jj=1,nElemsY
      DO ii=1,nElemsX


#ifdef MULTIFLUID
        !*If interface do not touch the conserved variables
        IF ( (Troubled_Cell_U(ii,jj) .NE. 0) ) THEN
          ! U(2:3,ii,jj)=0.5*(QtempM(2:3)+QtempP(2:3)) !*DO IT ONLY ON VELOCITY
          CYCLE
        END IF
#endif

        U(1:nVar,ii,jj)=0.5*(Cons_W_X(1:nVar,ii  ,jj)+Cons_W_X(1:nVar,ii+1,jj))+MESH_DX(1)/8.0*(slopes_W_X(1:nVar,ii  ,jj)-slopes_W_X(1:nVar,ii+1,jj))

      END DO
    END DO

  CASE DEFAULT
    ErrorMessage = "AF_PostProcessing_Potentially_HO_X not implemented for such reconstruction (order higher than 2)"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE AF_PostProcessing_Potentially_HO_X
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE AF_PostProcessing_Potentially_HO_Y()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_Equation,           ONLY: PrimToCons
USE MOD_Equation,           ONLY: ConsToPrim
USE MOD_Reconstruction,     ONLY: AF_MUSCL
USE MOD_Reconstruction,     ONLY: MINMOD_3_INPUTS
USE MOD_Reconstruction,     ONLY: SBMMUSCL
#ifdef MULTIFLUID
USE MOD_FiniteVolume2D_vars,ONLY: Which_Fluid_In_Cell_U
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_U
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_W_X
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_W_Y
USE MOD_Reconstruction,     ONLY: SIDED_AF_MUSCL

#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj, iVar
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: QtempM(1:nVar)
REAL               :: QtempP(1:nVar)
REAL               :: Qtemp(1:nVar), Vtemp(1:nVar)
REAL               :: Qaverage(1:nVar)
INTEGER            :: final_index
INTEGER            :: fluid_C, fluid_M, fluid_P
CHARACTER          :: which_side
REAL               :: slopes_U(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
REAL               :: slopes_W_Y( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
REAL               :: Cons_W_Y( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
REAL               :: Q_vector_dummy( 1:nVar,-1:+1)
REAL               :: theta=1.3
REAL               :: tau=0.5
REAL               :: WM_dummy(1:nVar,1:nGPs)
REAL               :: WP_dummy(1:nVar,1:nGPs)
!-------------------------------------------------------------------------------!

#ifdef MULTIFLUID

#ifdef LEVELSETINU
  final_index=nVar
#else
  final_index=nVar-1
#endif

#else

final_index=nVar

#endif


WM=0.0
WP=0.0

SELECT CASE (ABS(Reconstruction))
  CASE(1,10,2,20,21,22,23,24,25,26,27,28,29)

    !*Reconstruction of the limited boundary values

    DO ii=1,nElemsX
      DO jj=-1,nElemsY+2 

        !*NB: Limiting taking place in Conserved variables
        CALL PrimToCons(W_Y(1:nVar,ii,jj  ),QtempM(1:nVar))
        CALL PrimToCons(W_Y(1:nVar,ii,jj+1),QtempP(1:nVar))


#ifdef SIDEDPOSTPROCESSING
#ifdef MULTIFLUID
        fluid_C=Which_Fluid_In_Cell_U(ii,jj)
        fluid_M=Which_Fluid_In_Cell_U(ii,jj-1)
        fluid_P=Which_Fluid_In_Cell_U(ii,jj+1)

        IF ((fluid_C .EQ. fluid_M) .AND. (fluid_C .EQ. fluid_P)) THEN
          which_side="C"
        ELSE IF ( (fluid_C .EQ. fluid_M) .AND. (fluid_C .NE. fluid_P) ) THEN
          which_side="M"
        ELSE IF ( (fluid_C .NE. fluid_M) .AND. (fluid_C .EQ. fluid_P) ) THEN
          which_side="P"
        ELSE
          which_side="R"        
        END IF

        ! CALL SIDED_AF_MUSCL(&
        !           U(1:nVar  ,ii  ,jj-1:jj+1), &
        !           QtempM(1:nVar),&
        !           QtempP(1:nVar),&
        !           WM(1:nVar,1:nGPs,ii,jj),&
        !           WP(1:nVar,1:nGPs,ii,jj),&
        !           MESH_DX(2),which_side)

        PRINT*, "Sided PP not coded yet for this PP"
        STOP       
#else
        PRINT*, "Sided PP not existing without multilfuid flag"
        STOP
#endif
#else
        Q_vector_dummy(1:nVar,-1)=QtempM(1:nVar)
        Q_vector_dummy(1:nVar, 0)=U(1:nVar,ii,jj)
        Q_vector_dummy(1:nVar,+1)=QtempP(1:nVar)


#ifdef OVERCOMPRESSIVEPOSTPROCESSING
        CALL SBMMUSCL(Q_vector_dummy(1:nVar,-1:+1),WM_dummy(1:nVar,1:nGPs),WP_dummy(1:nVar,1:nGPs),MESH_DX(2)/2.0,slopes_U(1:nVar,ii,jj),tau_input=tau,theta_input=theta)
#else
        DO iVar=1,nVar
          slopes_U(iVar,ii,jj)=MINMOD_3_INPUTS((Q_vector_dummy(iVar, 0)-Q_vector_dummy(iVar,-1))/(0.5*MESH_DX(2))*theta, &
                                            &  (Q_vector_dummy(iVar,+1)-Q_vector_dummy(iVar,-1))/(    MESH_DX(2)), &
                                            &  (Q_vector_dummy(iVar,+1)-Q_vector_dummy(iVar, 0))/(0.5*MESH_DX(2))*theta &
                                            &  )
        END DO
#endif

#endif

      END DO
    END DO

    !*Replacement of the boundary values

    DO ii=1,nElemsX
      DO jj=0,nElemsY+2 

        Cons_W_Y(1:nVar,ii,jj)=0.5*(U(1:nVar,ii,jj-1)+U(1:nVar,ii  ,jj))+MESH_DX(1)/8.0*(slopes_U(1:nVar,ii,jj-1)-slopes_U(1:nVar,ii  ,jj))
        !*from the previous cell
        !*from the actual   cell

#ifdef MULTIFLUID
#ifdef SIDEDPOSTPROCESSING
        IF (Which_Fluid_In_Cell_U(ii,jj-1) .NE. Which_Fluid_In_Cell_U(ii,jj)) THEN
          CYCLE
        END IF
#else
        !*If interface do not touch the primitive variables
        IF (Troubled_Cell_W_Y(ii,jj) .NE. 0) THEN
          !*-----------------------
          !*SAFETY CHECK
          !*-----------------------
          ! IF ( (Troubled_Cell_U(ii,jj-1) .EQ. 0) .AND. (Troubled_Cell_U(ii,jj) .EQ. 0) ) THEN
          !   PRINT*, ii, jj
          !   PRINT*, Troubled_Cell_W_Y(ii,jj)
          !   PRINT*, Troubled_Cell_U(ii,jj-1), Troubled_Cell_U(ii,jj)
          !   STOP
          ! END IF 
          !*-----------------------

          CYCLE
        END IF
#endif
        !*-----------------------
        !*SAFETY CHECK
        !*-----------------------
        ! IF (Troubled_Cell_W_Y(ii,jj) .NE. 1) THEN
        !   IF ( (Troubled_Cell_U(ii,jj-1) .NE. 0) .OR. (Troubled_Cell_U(ii,jj) .NE. 0) ) THEN
        !     PRINT*, ii, jj
        !     PRINT*, Troubled_Cell_W_Y(ii,jj)
        !     PRINT*, Troubled_Cell_U(ii,jj-1), Troubled_Cell_U(ii,jj)
        !     STOP
        !   END IF 
        ! END IF 
        !*-----------------------
#endif

        CALL ConsToPrim(Cons_W_Y(1:nVar,ii,jj), W_Y(1:nVar,ii,jj)) 

      END DO
    END DO

    !*Slopes of the staggered mesh
    DO ii=1,nElemsX
      DO jj=1,nElemsY+1

        Q_vector_dummy(1:nVar,-1)=Cons_W_Y(1:nVar,ii,jj-1)
        Q_vector_dummy(1:nVar, 0)=Cons_W_Y(1:nVar,ii,jj)
        Q_vector_dummy(1:nVar,+1)=Cons_W_Y(1:nVar,ii,jj+1)

#ifdef OVERCOMPRESSIVEPOSTPROCESSING
        CALL SBMMUSCL(Q_vector_dummy(1:nVar,-1:+1),WM_dummy(1:nVar,1:nGPs),WP_dummy(1:nVar,1:nGPs),MESH_DX(2),slopes_W_Y(1:nVar,ii,jj),tau_input=tau,theta_input=theta)
#else
        DO iVar=1,nVar
          slopes_W_Y(iVar,ii,jj)=MINMOD_3_INPUTS((Q_vector_dummy(iVar, 0)-Q_vector_dummy(iVar,-1))/(    MESH_DX(2))*theta, &
                                              &  (Q_vector_dummy(iVar,+1)-Q_vector_dummy(iVar,-1))/(2.0*MESH_DX(2)), &
                                              &  (Q_vector_dummy(iVar,+1)-Q_vector_dummy(iVar, 0))/(    MESH_DX(2))*theta &
                                              &  )

        END DO
#endif

      END DO
    END DO




    !*Replacement of the conserved values
    DO ii=1,nElemsX
      DO jj=1,nElemsY

        CALL PrimToCons(W_Y(1:nVar,ii,jj  ),QtempM(1:nVar))
        CALL PrimToCons(W_Y(1:nVar,ii,jj+1),QtempP(1:nVar))

#ifdef MULTIFLUID
        !*If interface do not touch the conserved variables
        IF ( (Troubled_Cell_U(ii,jj) .NE. 0) ) THEN
          ! U(2:3,ii,jj)=0.5*(QtempM(2:3)+QtempP(2:3)) !*DO IT ONLY ON VELOCITY
          CYCLE
        END IF
#endif

        U(1:nVar,ii,jj)=0.5*(Cons_W_Y(1:nVar,ii  ,jj)+Cons_W_Y(1:nVar,ii,jj+1))+MESH_DX(2)/8.0*(slopes_W_Y(1:nVar,ii  ,jj)-slopes_W_Y(1:nVar,ii,jj+1))

      END DO
    END DO

  CASE DEFAULT
    ErrorMessage = "AF_PostProcessing_Potentially_HO_Y not implemented for such reconstruction (order higher than 2)"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE AF_PostProcessing_Potentially_HO_Y
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE Divergence_Free_Preparation_Of_Data_W_X()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX

USE MOD_FiniteVolume2D_vars,ONLY: W_X

USE MOD_FiniteVolume2D_vars,ONLY: NRows_X
USE MOD_FiniteVolume2D_vars,ONLY: From_2Indices_To_1Index_X
USE MOD_FiniteVolume2D_vars,ONLY: From_1Index_To_2Indices_X
USE MOD_FiniteVolume2D_vars,ONLY: Neighbours_1Index_X


USE MOD_Reconstruction       ,ONLY: First_Derivative_Central_Order2
USE MOD_JacobiIteration      ,ONLY: jacobi_2
USE MOD_IterativeLinearSolver,ONLY: pmgmres_ilu_cr
USE MOD_IterativeLinearSolver,ONLY: mgmres_st
use MOD_lsqr                 ,ONLY: lsqr_solver_ez
use algebra                  ,ONLY: Inverse
use MOD_DirectLinearSolver   ,ONLY: r8mat_fs
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER                             :: NNZsparse
INTEGER                             :: NRows
REAL   ,ALLOCATABLE, DIMENSION(:)   :: Values
INTEGER,ALLOCATABLE, DIMENSION(:)   :: Columns
INTEGER,ALLOCATABLE, DIMENSION(:)   :: Rows
INTEGER,ALLOCATABLE, DIMENSION(:)   :: RowStart
REAL   ,ALLOCATABLE, DIMENSION(:)   :: Diagonal
REAL   ,ALLOCATABLE, DIMENSION(:)   :: rhs
REAL   ,ALLOCATABLE, DIMENSION(:)   :: sol
REAL   ,ALLOCATABLE, DIMENSION(:)   :: sol_debug
REAL   ,ALLOCATABLE, DIMENSION(:)   :: sol_debug_2
REAL   ,ALLOCATABLE, DIMENSION(:)   :: sol_debug_3
REAL   ,ALLOCATABLE, DIMENSION(:)   :: sol_debug_4
REAL   ,ALLOCATABLE, DIMENSION(:,:) :: Mat
REAL   ,ALLOCATABLE, DIMENSION(:,:) :: invMat
REAL   ,ALLOCATABLE, DIMENSION(:)   :: initial_guess
REAL   ,ALLOCATABLE, DIMENSION(:,:) :: psi

INTEGER         :: ii, jj, iterations, its, indi
INTEGER         :: inde, indc, ind, indB, indR, indU, indL, indrow
REAL            :: dx, dy, Const_dx, Const_dy
REAL            :: dxvx, dyvy, dxv, dyu
REAL            :: tolerance=1e-14

type(lsqr_solver_ez) :: solver  !! main solver class
integer :: istop  !! solver exit code
!-------------------------------------------------------------------------------!


!*------------------------------
!*More or less, what we do here has been done for the IMEX system
!*In the end, we want to solve a laplacian (potential) problem
!*------------------------------

dx = MESH_DX(1)
dy = MESH_DX(2)

!*We already have NRows_X=(nElemsX+1)*nElemsY+2*(nElemsX+1)+2*nElemsY

!*NB: This is strongly low order because a 5 point stencil is assumed
!*5=C,D,R,U,L
NNZsparse = (nElemsX+1)*nElemsY*5 !*NB:+1 in X direction 

!*We need to take BCs into account now
!*We assume strong
!------------------------------!
! Bottom Boundary Conditions   !
!------------------------------!
NNZsparse = NNZsparse + nElemsX+1

!------------------------------!
! Right Boundary Conditions    !
!------------------------------!
NNZsparse = NNZsparse + nElemsY

!------------------------------!
! Top Boundary Conditions      !
!------------------------------!
NNZsparse = NNZsparse + nElemsX+1

!------------------------------!
! Left Boundary Conditions     !
!------------------------------!
NNZsparse = NNZsparse + nElemsY

ALLOCATE(Values(NNZsparse)                                            )
ALLOCATE(Columns(NNZsparse)                                           )
ALLOCATE(Rows(NNZsparse)                                              )
ALLOCATE(RowStart(NRows_X+1)                                          )
ALLOCATE(Diagonal(NRows_X)                                            )
ALLOCATE(rhs(NRows_X)                                                 )
ALLOCATE(sol(NRows_X)                                                 )
ALLOCATE(initial_guess(NRows_X)                                       )
ALLOCATE(psi(-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1) )


Values     = 0.0
Columns    = 0
Rows       = 0
RowStart   = 0
Diagonal   = 0.0
rhs        = 0.0
sol        = 0.0




!*===========================
!*MATRIX
!*===========================

!*-------------------------------
!*Switching to rescaling by dx**2
!*-------------------------------
#if defined(NORMALIZELINEARSYSTEMS) || defined(NORMALIZEALL)
Const_dx = 1.0
Const_dy = (dx/dy)**2
#else
Const_dx = 1.0/dx**2
Const_dy = 1.0/dy**2
#endif
!*-------------------------------


!*First equations->Inside the domain

indc=0 !*index cell
inde=0 !*index element
indrow=0 !*index row

DO jj=1,nElemsY
  DO ii=1,nElemsX+1 

    !*C,B,R,U,L
    indc=indc+1
    inde=inde+1
    indrow=indrow+1

    RowStart(indc)=inde

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
    Columns(inde)  = indc
    Rows(inde)     = indrow
    Values(inde)   = -2.0*Const_dx-2.0*Const_dy
    Diagonal(indc) = Values(inde)

    !*B (ii,jj-1)
    inde=inde+1
    indB=Neighbours_1Index_X(indc,2) !*From_2Indices_To_1Index_X(ii,jj-1)
    Columns(inde) = indB
    Rows(inde)     = indrow
    Values(inde)  = +Const_dy

    !*R (ii+1,jj)
    inde=inde+1
    indR=Neighbours_1Index_X(indc,3) !*From_2Indices_To_1Index_X(ii+1,jj)
    Columns(inde) = indR
    Rows(inde)     = indrow
    Values(inde)  = +Const_dx

    !*U (ii,jj+1)
    inde=inde+1
    indU=Neighbours_1Index_X(indc,4) !*From_2Indices_To_1Index_X(ii,jj+1)
    Columns(inde) = indU
    Rows(inde)     = indrow
    Values(inde)  = +Const_dy

    !*L (ii-1,jj)
    inde=inde+1
    indL=Neighbours_1Index_X(indc,5) !*From_2Indices_To_1Index_X(ii-1,jj)
    Columns(inde) = indL
    Rows(inde)     = indrow
    Values(inde)  = +Const_dx

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
  indrow=indrow+1
  RowStart(indc) = inde
  Columns(inde)  = indc
  Rows(inde)  = indrow
  Values(inde)   = 1.0
  Diagonal(indc) = Values(inde)
  !*Nothing u=something
END DO


!*->Right boundary
ii=nElemsX+2
DO jj=1,nElemsY
  indc=indc+1
  inde=inde+1
  indrow=indrow+1
  RowStart(indc) = inde
  Columns(inde)  = indc
  Rows(inde)  = indrow
  Values(inde)   = 1.0
  Diagonal(indc) = Values(inde)
  !*Nothing u=something
END DO

!*->Top boundary
jj=nElemsY+1
DO ii=1,nElemsX+1 
  indc=indc+1
  inde=inde+1
  indrow=indrow+1
  RowStart(indc) = inde
  Columns(inde)  = indc
  Rows(inde)  = indrow
  Values(inde)   = 1.0
  Diagonal(indc) = Values(inde)
  !*Nothing u=something
END DO



!*->Left boundary
ii=0
DO jj=1,nElemsY
  indc=indc+1
  inde=inde+1
  indrow=indrow+1
  RowStart(indc) = inde
  Columns(inde)  = indc
  Rows(inde)  = indrow
  Values(inde)   = 1.0
  Diagonal(indc) = Values(inde)
  !*Nothing u=something
END DO


RowStart(NRows_X+1)=NNZsparse+1

!*-----------------------------------------
!*SAFETY CHECK
!*-----------------------------------------
! PRINT*, Columns
! PRINT*, Values
! PRINT*, RowStart
! PRINT*, Diagonal
! STOP
!*-----------------------------------------




!*===========================
!*RIGHT-HAND SIDE
!*===========================
rhs=0.0

!*Let us start by the known contribution inside the domain
!*Inside the domain
indc=0
DO jj=1,nElemsY
  DO ii=1,nElemsX+1 !*NB: +1 in X direction
    indc=indc+1

    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc

    !*NB: It is important to have BCs here
    ! dxvx=First_Derivative_Central_Order2(W_X(2,ii-1:ii+1,jj),dx)
    ! dyvy=First_Derivative_Central_Order2(W_X(3,ii,jj-1:jj+1),dy)

    dxv=First_Derivative_Central_Order2(W_X(3,ii-1:ii+1,jj),dx)
    dyu=First_Derivative_Central_Order2(W_X(2,ii,jj-1:jj+1),dy)

    !*-------------------------------
    !*Switching to rescaling by dx**2
    !*-------------------------------
#if defined(NORMALIZELINEARSYSTEMS) || defined(NORMALIZEALL)
    rhs(indc)=(-dxv+dyu)*dx**2
#else
    rhs(indc)=-dxv+dyu
#endif
    !*-------------------------------
                          
  END DO !ii
END DO !jj



!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX+1 
  indc=indc+1
  !*Do nothing, keep 0
END DO


!*->Right boundary
ii=nElemsX+2
DO jj=1,nElemsY
  indc=indc+1
  !*Do nothing, keep 0
END DO

!*->Top boundary
jj=nElemsY+1
DO ii=1,nElemsX+1 
  indc=indc+1
  !*Do nothing, keep 0
END DO



!*->Left boundary
ii=0
DO jj=1,nElemsY
  indc=indc+1
  !*Do nothing, keep 0
END DO


!*===========================
!*SOLUTION
!*===========================
its = 10000 
sol = 0.0
call pmgmres_ilu_cr(NRows_X, NNZsparse, RowStart, Columns, Values, sol, rhs, its, 50, 1e-12, 1e-8 )


#if(1==0)
PRINT*, "TEST 1 Jacobi"
ALLOCATE(sol_debug(NRows_X))
initial_guess=0.0
iterations=0
CALL jacobi_2(RowStart,Columns,Values,Diagonal,rhs,sol_debug,iterations,initial_guess)
DO indi=1,NRows_X
  IF(ABS(sol(indi)-sol_debug(indi)) .GT. tolerance) THEN
    PRINT*, "Jacobi and the other solver give different results in W_X"
    PRINT*, indi, sol(indi)-sol_debug(indi)
    STOP
  END IF
END DO
DEALLOCATE(sol_debug     )
#endif

#if(1==0)
PRINT*, "TEST 2 other iterative solver"
ALLOCATE(sol_debug_2(NRows_X))
CALL solver%initialize(NRows_X,NRows_X,Values,Rows,Columns) ! use defaults for other optional inputs
CALL solver%solve(rhs,0.0,sol_debug_2,istop)       ! solve the linear system
WRITE(*,*) 'istop = ', istop
DO indi=1,NRows_X
  IF(ABS(sol(indi)-sol_debug_2(indi)) .GT. tolerance) THEN
    PRINT*, "New solver and the other solver give different results in W_X"
    PRINT*, indi, sol(indi)-sol_debug_2(indi)
  END IF
END DO
DEALLOCATE(sol_debug_2   )
#endif


#if(1==0)
PRINT*, "TEST 3 direct solver"
ALLOCATE(sol_debug_3(NRows_X))
ALLOCATE(Mat(NRows_X,NRows_X))
Mat=0.0
DO inde=1,NNZsparse
  Mat(Rows(inde),Columns(inde))=Values(inde)
END DO
CALL r8mat_fs ( NRows_X, Mat, rhs, sol_debug_3 )
DO indi=1,NRows_X
  IF(ABS(sol(indi)-sol_debug_3(indi)) .GT. tolerance) THEN
    PRINT*, "New solver and direct solver give different results in W_X"
    PRINT*, indi, sol(indi)-sol_debug_3(indi)
    STOP
  END IF
END DO
DEALLOCATE(sol_debug_3   )
DEALLOCATE(Mat           )
#endif


#if(1==0)
PRINT*, "TEST 4 GMRES"
ALLOCATE(sol_debug_4(NRows_X))
its = 10000 
sol_debug_4 = 0.0
CALL mgmres_st(NRows_X, NNZsparse, Rows, Columns, Values, sol_debug_4, rhs, its, 5, 1e-14, 1e-14 )
DO indi=1,NRows_X
  IF(ABS(sol(indi)-sol_debug_4(indi)) .GT. tolerance) THEN
    PRINT*, "GMRES and the other solver give different results in W_X"
    PRINT*, indi, sol(indi)-sol_debug_4(indi)
    STOP
  END IF
END DO
DEALLOCATE(sol_debug_4   )
STOP
#endif

!*===========================
!*PUT IN MATRIX
!*===========================
indc=0
!*Inside domain
DO jj=1,nElemsY
  DO ii=1,nElemsX+1 !*NB:+1 in X direction
    indc=indc+1
    psi(ii,jj)=sol(indc)
  END DO
END DO

!*Ghosts
!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX+1 
  indc=indc+1
  psi(ii,jj)=sol(indc)
END DO

!*->Right boundary
ii=nElemsX+2
DO jj=1,nElemsY
  indc=indc+1
  psi(ii,jj)=sol(indc)
END DO

!*->Top boundary
jj=nElemsY+1
DO ii=1,nElemsX+1 
  indc=indc+1
  psi(ii,jj)=sol(indc)
END DO

!*->Left boundary
ii=0
DO jj=1,nElemsY
  indc=indc+1
  psi(ii,jj)=sol(indc)
END DO

!*===========================
!*OVERWRITE W_X IN VELOCITY FIELD
!*===========================
DO jj=1,nElemsY
  DO ii=1,nElemsX+1 !*NB:+1 in X direction
    ! W_X(2,ii,jj)=First_Derivative_Central_Order2(psi(ii-1:ii+1,jj),dx)
    ! W_X(3,ii,jj)=-First_Derivative_Central_Order2(psi(ii,jj-1:jj+1),dy)
    W_X(2,ii,jj)=First_Derivative_Central_Order2(psi(ii,jj-1:jj+1),dy)
    W_X(3,ii,jj)=-First_Derivative_Central_Order2(psi(ii-1:ii+1,jj),dx)
  END DO
END DO


! DO jj=1,nElemsY
!   DO ii=1,nElemsX+1 !*NB:+1 in X direction
!     PRINT*, ii, jj, First_Derivative_Central_Order2(W_X(2,ii-1:ii+1,jj),dx)+First_Derivative_Central_Order2(W_X(3,ii,jj-1:jj+1),dy)
!   END DO
! END DO



DEALLOCATE(Values        )
DEALLOCATE(Columns       )
DEALLOCATE(Rows          )
DEALLOCATE(RowStart      )
DEALLOCATE(Diagonal      )
DEALLOCATE(rhs           )
DEALLOCATE(sol           )
DEALLOCATE(initial_guess )
DEALLOCATE(psi           )




END SUBROUTINE Divergence_Free_Preparation_Of_Data_W_X
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE Divergence_Free_Preparation_Of_Data_W_Y()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX

USE MOD_FiniteVolume2D_vars,ONLY: W_Y

USE MOD_FiniteVolume2D_vars,ONLY: NRows_Y
USE MOD_FiniteVolume2D_vars,ONLY: From_2Indices_To_1Index_Y
USE MOD_FiniteVolume2D_vars,ONLY: From_1Index_To_2Indices_Y
USE MOD_FiniteVolume2D_vars,ONLY: Neighbours_1Index_Y


USE MOD_Reconstruction     ,ONLY: First_Derivative_Central_Order2
USE MOD_JacobiIteration    ,ONLY: jacobi_2
USE MOD_IterativeLinearSolver, ONLY: pmgmres_ilu_cr
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER                             :: NNZsparse
INTEGER                             :: NRows
REAL   ,ALLOCATABLE, DIMENSION(:)   :: Values
INTEGER,ALLOCATABLE, DIMENSION(:)   :: Columns
INTEGER,ALLOCATABLE, DIMENSION(:)   :: RowStart
REAL   ,ALLOCATABLE, DIMENSION(:)   :: Diagonal
REAL   ,ALLOCATABLE, DIMENSION(:)   :: rhs
REAL   ,ALLOCATABLE, DIMENSION(:)   :: sol
REAL   ,ALLOCATABLE, DIMENSION(:)   :: sol_debug
REAL   ,ALLOCATABLE, DIMENSION(:)   :: initial_guess
REAL   ,ALLOCATABLE, DIMENSION(:,:) :: psi

INTEGER         :: ii, jj, iterations, its, indi
INTEGER         :: inde, indc, ind, indB, indR, indU, indL
REAL            :: dx, dy, Const_dx, Const_dy
REAL            :: dxvx, dyvy, dxv, dyu
REAL            :: tolerance=1e-9
!-------------------------------------------------------------------------------!


!*------------------------------
!*More or less, what we do here has been done for the IMEX system
!*In the end, we want to solve a laplacian (potential) problem
!*------------------------------

dx = MESH_DX(1)
dy = MESH_DX(2)

!*We already have NRows_Y=nElemsX*(nElemsY+1)+2*nElemsX+2*(nElemsY+1)

!*NB: This is strongly low order because a 5 point stencil is assumed
!*5=C,D,R,U,L
NNZsparse = nElemsX*(nElemsY+1)*5 !*NB:+1 in Y direction 

!*We need to take BCs into account now
!*We assume strong
!------------------------------!
! Bottom Boundary Conditions   !
!------------------------------!
NNZsparse = NNZsparse + nElemsX

!------------------------------!
! Right Boundary Conditions    !
!------------------------------!
NNZsparse = NNZsparse + nElemsY+1

!------------------------------!
! Top Boundary Conditions      !
!------------------------------!
NNZsparse = NNZsparse + nElemsX

!------------------------------!
! Left Boundary Conditions     !
!------------------------------!
NNZsparse = NNZsparse + nElemsY+1

ALLOCATE(Values(NNZsparse)                                            )
ALLOCATE(Columns(NNZsparse)                                           )
ALLOCATE(RowStart(NRows_Y+1)                                          )
ALLOCATE(Diagonal(NRows_Y)                                            )
ALLOCATE(rhs(NRows_Y)                                                 )
ALLOCATE(sol(NRows_Y)                                                 )
ALLOCATE(sol_debug(NRows_Y)                                                 )
ALLOCATE(initial_guess(NRows_Y)                                       )
ALLOCATE(psi(-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1) )


Values     = 0.0
Columns    = 0
RowStart   = 0
Diagonal   = 0.0
rhs        = 0.0
sol        = 0.0




!*===========================
!*MATRIX
!*===========================


!*-------------------------------
!*Switching to rescaling by dx**2
!*-------------------------------
#if defined(NORMALIZELINEARSYSTEMS) || defined(NORMALIZEALL)
Const_dx = 1.0
Const_dy = (dx/dy)**2
#else
Const_dx = 1.0/dx**2
Const_dy = 1.0/dy**2
#endif
!*-------------------------------

!*First equations->Inside the domain

indc=0 !*index cell
inde=0 !*index element

DO jj=1,nElemsY+1
  DO ii=1,nElemsX

    !*C,B,R,U,L
    indc=indc+1
    inde=inde+1

    RowStart(indc)=inde

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
    ! ind=Neighbours_1Index_Y(indc,1) !*NB: This is equal to indc
    !*-----------------------------------------
    !*As a safety check this should give ii and jj
    !*-----------------------------------------
    ! IF ((indc-ind) .NE. 0) THEN
    !   PRINT*, "Problem"
    !   PRINT*, ind-indc
    !   STOP
    ! END IF
    !*-----------------------------------------
    Columns(inde)  = indc
    Values(inde)   = -2.0*Const_dx-2.0*Const_dy
    Diagonal(indc) = Values(inde)

    !*B (ii,jj-1)
    inde=inde+1
    indB=Neighbours_1Index_Y(indc,2) !*From_2Indices_To_1Index_X(ii,jj-1)
    Columns(inde) = indB
    Values(inde)  = +Const_dy

    !*R (ii+1,jj)
    inde=inde+1
    indR=Neighbours_1Index_Y(indc,3) !*From_2Indices_To_1Index_X(ii+1,jj)
    Columns(inde) = indR
    Values(inde)  = +Const_dx

    !*U (ii,jj+1)
    inde=inde+1
    indU=Neighbours_1Index_Y(indc,4) !*From_2Indices_To_1Index_X(ii,jj+1)
    Columns(inde) = indU
    Values(inde)  = +Const_dy

    !*L (ii-1,jj)
    inde=inde+1
    indL=Neighbours_1Index_Y(indc,5) !*From_2Indices_To_1Index_X(ii-1,jj)
    Columns(inde) = indL
    Values(inde)  = +Const_dx

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
  RowStart(indc) = inde
  Columns(inde)  = indc
  Values(inde)   = 1.0
  Diagonal(indc) = Values(inde)
  !*Nothing u=something
END DO


!*->Right boundary
ii=nElemsX+1
DO jj=1,nElemsY+1
  indc=indc+1
  inde=inde+1
  RowStart(indc) = inde
  Columns(inde)  = indc
  Values(inde)   = 1.0
  Diagonal(indc) = Values(inde)
  !*Nothing u=something
END DO

!*->Top boundary
jj=nElemsY+2
DO ii=1,nElemsX 
  indc=indc+1
  inde=inde+1
  RowStart(indc) = inde
  Columns(inde)  = indc
  Values(inde)   = 1.0
  Diagonal(indc) = Values(inde)
  !*Nothing u=something
END DO



!*->Left boundary
ii=0
DO jj=1,nElemsY+1
  indc=indc+1
  inde=inde+1
  RowStart(indc) = inde
  Columns(inde)  = indc
  Values(inde)   = 1.0
  Diagonal(indc) = Values(inde)
  !*Nothing u=something
END DO


RowStart(NRows_Y+1)=NNZsparse+1

!*-----------------------------------------
!*SAFETY CHECK
!*-----------------------------------------
! PRINT*, Columns
! PRINT*, Values
! PRINT*, RowStart
! PRINT*, Diagonal
! STOP
!*-----------------------------------------




!*===========================
!*RIGHT-HAND SIDE
!*===========================
rhs=0.0

!*Let us start by the known contribution inside the domain
!*Inside the domain
indc=0
DO jj=1,nElemsY+1 !*NB: +1 in Y direction
  DO ii=1,nElemsX
    indc=indc+1

    !*NB: We are sweeping the elements in the same order we gave
    !*Thus, we do not need to retrive the global number of ii,jj.
    !*It is indc

    !*NB: It is important to have BCs here
    ! dxvx=First_Derivative_Central_Order2(W_Y(2,ii-1:ii+1,jj),dx)
    ! dyvy=First_Derivative_Central_Order2(W_Y(3,ii,jj-1:jj+1),dy)

    dxv=First_Derivative_Central_Order2(W_Y(3,ii-1:ii+1,jj),dx)
    dyu=First_Derivative_Central_Order2(W_Y(2,ii,jj-1:jj+1),dy)

    !*-------------------------------
    !*Switching to rescaling by dx**2
    !*-------------------------------
#if defined(NORMALIZELINEARSYSTEMS) || defined(NORMALIZEALL)
    rhs(indc)=(-dxv+dyu)*dx**2
#else
    rhs(indc)=-dxv+dyu
#endif
    !*-------------------------------


  END DO !ii
END DO !jj



!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX 
  indc=indc+1
  !*Do nothing, keep 0
END DO


!*->Right boundary
ii=nElemsX+1
DO jj=1,nElemsY+1
  indc=indc+1
  !*Do nothing, keep 0
END DO

!*->Top boundary
jj=nElemsY+2
DO ii=1,nElemsX 
  indc=indc+1
  !*Do nothing, keep 0
END DO



!*->Left boundary
ii=0
DO jj=1,nElemsY+1
  indc=indc+1
  !*Do nothing, keep 0
END DO


!*===========================
!*SOLUTION
!*===========================
its = 10000 
sol = 0.0
call pmgmres_ilu_cr(NRows_Y, NNZsparse, RowStart, Columns, Values, sol, rhs, its, 50, 1e-12, 1e-8 )


#if(1==0)
initial_guess=0.0
iterations=0
CALL jacobi_2(RowStart,Columns,Values,Diagonal,rhs,sol_debug,iterations,initial_guess)
DO indi=1,NRows_Y
  IF(ABS(sol(indi)-sol_debug(indi)) .GT. tolerance) THEN
    PRINT*, "Jacobi and the other solver give different results in W_Y"
    PRINT*, indi, sol(indi)-sol_debug(indi)
    STOP
  END IF
END DO
#endif


!*===========================
!*PUT IN MATRIX
!*===========================
indc=0
!*Inside domain
DO jj=1,nElemsY+1 !*NB:+1 in Y direction
  DO ii=1,nElemsX
    indc=indc+1
    psi(ii,jj)=sol(indc)
  END DO
END DO

!*Ghosts
!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX 
  indc=indc+1
  psi(ii,jj)=sol(indc)
END DO

!*->Right boundary
ii=nElemsX+1
DO jj=1,nElemsY+1
  indc=indc+1
  psi(ii,jj)=sol(indc)
END DO

!*->Top boundary
jj=nElemsY+2
DO ii=1,nElemsX 
  indc=indc+1
  psi(ii,jj)=sol(indc)
END DO

!*->Left boundary
ii=0
DO jj=1,nElemsY+1
  indc=indc+1
  psi(ii,jj)=sol(indc)
END DO

!*===========================
!*OVERWRITE W_Y IN VELOCITY FIELD
!*===========================
DO jj=1,nElemsY+1 !*NB:+1 in Y direction
  DO ii=1,nElemsX
    ! W_Y(2,ii,jj)=First_Derivative_Central_Order2(psi(ii-1:ii+1,jj),dx)
    ! W_Y(3,ii,jj)=-First_Derivative_Central_Order2(psi(ii,jj-1:jj+1),dy)
    W_Y(2,ii,jj)=First_Derivative_Central_Order2(psi(ii,jj-1:jj+1),dy)
    W_Y(3,ii,jj)=-First_Derivative_Central_Order2(psi(ii-1:ii+1,jj),dx)
  END DO
END DO


! DO jj=1,nElemsY+1 !*NB:+1 in Y direction
!   DO ii=1,nElemsX
!     PRINT*, ii, jj, First_Derivative_Central_Order2(W_Y(2,ii-1:ii+1,jj),dx)+First_Derivative_Central_Order2(W_Y(3,ii,jj-1:jj+1),dy)
!   END DO
! END DO

DEALLOCATE(Values        )
DEALLOCATE(Columns       )
DEALLOCATE(RowStart      )
DEALLOCATE(Diagonal      )
DEALLOCATE(rhs           )
DEALLOCATE(sol           )
DEALLOCATE(sol_debug     )
DEALLOCATE(initial_guess )
DEALLOCATE(psi           )




END SUBROUTINE Divergence_Free_Preparation_Of_Data_W_Y
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
#if defined(IMEX) || defined(IMEXMOMENTUM)
SUBROUTINE Put_Matrix_In_Vector_X(M,vec)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: NRows_X
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, INTENT(IN)    :: M(-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
REAL, INTENT(INOUT) :: vec(1:NRows_X)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, indc
!-------------------------------------------------------------------------------!

indc=0
!*Inside domain
DO jj=1,nElemsY
  DO ii=1,nElemsX+1 !*NB:+1 in X direction
    indc=indc+1
    vec(indc)=M(ii,jj)
  END DO
END DO

!*Ghosts
!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX+1 
  indc=indc+1
  vec(indc)=M(ii,jj)
END DO

!*->Right boundary
ii=nElemsX+2
DO jj=1,nElemsY
  indc=indc+1
  vec(indc)=M(ii,jj)
END DO

!*->Top boundary
jj=nElemsY+1
DO ii=1,nElemsX+1 
  indc=indc+1
  vec(indc)=M(ii,jj)
END DO

!*->Left boundary
ii=0
DO jj=1,nElemsY
  indc=indc+1
  vec(indc)=M(ii,jj)
END DO


END SUBROUTINE Put_Matrix_In_Vector_X
#endif
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
#if defined(IMEX) || defined(IMEXMOMENTUM)
SUBROUTINE Put_Matrix_In_Vector_Y(M,vec)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: NRows_Y
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, INTENT(IN)    :: M(-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
REAL, INTENT(INOUT) :: vec(1:NRows_Y)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, indc
!-------------------------------------------------------------------------------!

indc=0
!*Inside domain
DO jj=1,nElemsY+1 !*NB:+1 in Y direction
  DO ii=1,nElemsX
    indc=indc+1
    vec(indc)=M(ii,jj)
  END DO
END DO

!*Ghosts
!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX 
  indc=indc+1
  vec(indc)=M(ii,jj)
END DO

!*->Right boundary
ii=nElemsX+1
DO jj=1,nElemsY+1
  indc=indc+1
  vec(indc)=M(ii,jj)
END DO

!*->Top boundary
jj=nElemsY+2
DO ii=1,nElemsX 
  indc=indc+1
  vec(indc)=M(ii,jj)
END DO

!*->Left boundary
ii=0
DO jj=1,nElemsY+1
  indc=indc+1
  vec(indc)=M(ii,jj)
END DO


END SUBROUTINE Put_Matrix_In_Vector_Y
#endif
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef CENTEREDPRIMITIVE
#if defined(IMEX) || defined(IMEXMOMENTUM)
SUBROUTINE Put_Matrix_In_Vector(M,vec)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: NRows_WC
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, INTENT(IN)    :: M(-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
REAL, INTENT(INOUT) :: vec(1:NRows_WC)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, indc
!-------------------------------------------------------------------------------!

indc=0
!*Inside domain
DO jj=1,nElemsY
  DO ii=1,nElemsX
    indc=indc+1
    vec(indc)=M(ii,jj)
  END DO
END DO

!*Ghosts
!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX
  indc=indc+1
  vec(indc)=M(ii,jj)
END DO

!*->Right boundary
ii=nElemsX+1
DO jj=1,nElemsY
  indc=indc+1
  vec(indc)=M(ii,jj)
END DO

!*->Top boundary
jj=nElemsY+1
DO ii=1,nElemsX 
  indc=indc+1
  vec(indc)=M(ii,jj)
END DO

!*->Left boundary
ii=0
DO jj=1,nElemsY
  indc=indc+1
  vec(indc)=M(ii,jj)
END DO


END SUBROUTINE Put_Matrix_In_Vector
#endif
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
#if defined(IMEX) || defined(IMEXMOMENTUM)
SUBROUTINE Put_Vector_In_Matrix_X(vec,M)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: NRows_X
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, INTENT(IN)       :: vec(1:NRows_X)
REAL, INTENT(INOUT)    :: M(-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, indc
!-------------------------------------------------------------------------------!

indc=0
!*Inside domain
DO jj=1,nElemsY
  DO ii=1,nElemsX+1 !*NB:+1 in X direction
    indc=indc+1
    M(ii,jj)=vec(indc)
  END DO
END DO

!*------------------------------
!*Added because I need ghosts for px and py
!*------------------------------
!*Ghosts
!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX+1 
  indc=indc+1
  M(ii,jj)=vec(indc)
END DO

!*->Right boundary
ii=nElemsX+2
DO jj=1,nElemsY
  indc=indc+1
  M(ii,jj)=vec(indc)
END DO

!*->Top boundary
jj=nElemsY+1
DO ii=1,nElemsX+1 
  indc=indc+1
  M(ii,jj)=vec(indc)
END DO

!*->Left boundary
ii=0
DO jj=1,nElemsY
  indc=indc+1
  M(ii,jj)=vec(indc)
END DO

END SUBROUTINE Put_Vector_In_Matrix_X
#endif
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
#if defined(IMEX) || defined(IMEXMOMENTUM)
SUBROUTINE Put_Vector_In_Matrix_Y(vec,M)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: NRows_Y
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, INTENT(IN)       :: vec(1:NRows_Y)
REAL, INTENT(INOUT)    :: M(-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, indc
!-------------------------------------------------------------------------------!

indc=0
!*Inside domain
DO jj=1,nElemsY+1 !*NB:+1 in Y direction
  DO ii=1,nElemsX
    indc=indc+1
    M(ii,jj)=vec(indc)
  END DO
END DO

!*------------------------------
!*Added because I need ghosts for px and py
!*------------------------------
!*Ghosts
!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX 
  indc=indc+1
  M(ii,jj)=vec(indc)
END DO

!*->Right boundary
ii=nElemsX+1
DO jj=1,nElemsY+1
  indc=indc+1
  M(ii,jj)=vec(indc)
END DO

!*->Top boundary
jj=nElemsY+2
DO ii=1,nElemsX 
  indc=indc+1
  M(ii,jj)=vec(indc)
END DO

!*->Left boundary
ii=0
DO jj=1,nElemsY+1
  indc=indc+1
  M(ii,jj)=vec(indc)
END DO


END SUBROUTINE Put_Vector_In_Matrix_Y
#endif
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef CENTEREDPRIMITIVE
#if defined(IMEX) || defined(IMEXMOMENTUM)
SUBROUTINE Put_Vector_In_Matrix(vec,M)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: NRows_WC
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, INTENT(IN)       :: vec(1:NRows_WC)
REAL, INTENT(INOUT)    :: M(-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, indc
!-------------------------------------------------------------------------------!

indc=0
!*Inside domain
DO jj=1,nElemsY
  DO ii=1,nElemsX
    indc=indc+1
    M(ii,jj)=vec(indc)
  END DO
END DO

!*------------------------------
!*Added because I need ghosts for px and py
!*------------------------------
!*Ghosts
!*Later equations->Boundary conditions
!*->Down boundary
jj=0
DO ii=1,nElemsX
  indc=indc+1
  M(ii,jj)=vec(indc)
END DO

!*->Right boundary
ii=nElemsX+1
DO jj=1,nElemsY
  indc=indc+1
  M(ii,jj)=vec(indc)
END DO

!*->Top boundary
jj=nElemsY+1
DO ii=1,nElemsX
  indc=indc+1
  M(ii,jj)=vec(indc)
END DO

!*->Left boundary
ii=0
DO jj=1,nElemsY
  indc=indc+1
  M(ii,jj)=vec(indc)
END DO

END SUBROUTINE Put_Vector_In_Matrix
#endif
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE FVTimeDerivativeConservedOnly(t,Optional_Input)
!-------------------------------------------------------------------------------!
USE MOD_Equation,       ONLY: BoundaryConditions
USE MOD_Equation,       ONLY: SourceTerms
USE MOD_Equation,       ONLY: SourceTerm_Conserved_INPUT_CONSERVED
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
USE MOD_Equation,       ONLY: SourceTerm_Primitive_INPUT_PRIMITIVE
#endif
USE MOD_Reconstruction, ONLY: ReconstructionX
USE MOD_Reconstruction, ONLY: ReconstructionY
USE MOD_Reconstruction, ONLY: ReconstructionFixX
USE MOD_Reconstruction, ONLY: ReconstructionFixY
USE MOD_ShocksIndicator,ONLY: ShocksIndicatorX
USE MOD_ShocksIndicator,ONLY: ShocksIndicatorY
#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars, ONLY: WC
USE MOD_FiniteVolume2D_vars, ONLY: WCt
#ifdef PATHCONSERVATIVESHOCKDETECTION
USE MOD_FiniteVolume2D_vars, ONLY: Cell_To_Limit_PCSD
#endif
#endif
#ifdef ACTIVEFLUX
USE MOD_Reconstruction, ONLY: FromCellAveragesToPointValues

!*Unified subroutines
USE MOD_FiniteVolume2D_vars, ONLY: W_X
USE MOD_FiniteVolume2D_vars, ONLY: W_Y
USE MOD_FiniteVolume2D_vars, ONLY: Wt_X
USE MOD_FiniteVolume2D_vars, ONLY: Wt_Y
#endif
USE MOD_Reconstruction,      ONLY: ReconstructionX_INPUT_PRIMITIVE
USE MOD_Reconstruction,      ONLY: ReconstructionY_INPUT_PRIMITIVE
USE MOD_Reconstruction,      ONLY: ReconstructionX_INPUT_CONSERVED
USE MOD_Reconstruction,      ONLY: ReconstructionY_INPUT_CONSERVED
USE MOD_Reconstruction,      ONLY: ReconstructionFixX_INPUT_PRIMITIVE
USE MOD_Reconstruction,      ONLY: ReconstructionFixY_INPUT_PRIMITIVE
USE MOD_Reconstruction,      ONLY: ReconstructionFixX_INPUT_CONSERVED
USE MOD_Reconstruction,      ONLY: ReconstructionFixY_INPUT_CONSERVED
USE MOD_FiniteVolume2D_vars, ONLY: nGhosts
USE MOD_FiniteVolume2D_vars, ONLY: nElemsX
USE MOD_FiniteVolume2D_vars, ONLY: nElemsY
USE MOD_FiniteVolume2D_vars, ONLY: U
USE MOD_FiniteVolume2D_vars, ONLY: nVar
#ifdef MULTIFLUID
USE MOD_FiniteVolume2D_vars, ONLY: Troubled_Cell_W_X
USE MOD_FiniteVolume2D_vars, ONLY: Troubled_Cell_W_Y
#endif
USE MOD_FiniteVolume2D_vars, ONLY: nGhosts
USE MOD_FiniteVolume2D_vars, ONLY: nGPs
USE MOD_FiniteVolume2D_vars, ONLY: WM
USE MOD_FiniteVolume2D_vars, ONLY: WP
USE MOD_FiniteVolume2D_vars, ONLY: MeshBary
USE MOD_FiniteVolume2D_vars, ONLY: MeshGP
USE MOD_FiniteVolume2D_vars, ONLY: nDims
#ifdef ACTIVEFLUX
USE MOD_FiniteVolume2D_vars, ONLY: MeshBary_X
USE MOD_FiniteVolume2D_vars, ONLY: MeshGP_X
USE MOD_FiniteVolume2D_vars, ONLY: MeshBary_Y
USE MOD_FiniteVolume2D_vars, ONLY: MeshGP_Y
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,    INTENT(IN)           :: t
INTEGER, INTENT(IN), OPTIONAL :: Optional_Input
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
#ifdef CENTEREDPRIMITIVE
REAL            :: W_support( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) !*NB: Exactly as U
#endif
#ifdef ACTIVEFLUX
REAL            :: W_X_support( 1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)!*NB:+1 in X direction with respect to U
REAL            :: W_Y_support( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)!*NB:+1 in Y direction with respect to U
#endif
!*For passing from conserved to primitive
REAL            :: WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1) 
REAL            :: WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1) 
INTEGER         :: ii,jj


CALL BoundaryConditions(t)




#if defined(CENTEREDPRIMITIVE)
!*=========================
!*U
!*=========================
#ifdef DETECTEACHSTAGE
IF (.NOT.(PRESENT(Optional_Input))) THEN
  CALL Shock_Detector_Based_On_Path_Conservative()
  ! PRINT*, "Detect"
ELSE
  ! PRINT*, "Not detect"
END IF
#endif


#ifdef RECONSTRUCTFROMCONSERVED
CALL ComputeTransitionMatricesConservedInput(U( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),  &  !*Vector
                          & -nGhosts,nElemsX+nGhosts+1, &  !*X bounds vector  
                          & -nGhosts,nElemsY+nGhosts+1,   &  !*Y bounds vector
                          & 1,nElemsX, &                   !*X bounds loop    
                          & 1,nElemsY)                       !*Y bounds loop 
#else
CALL ComputeTransitionMatricesPrimitiveInput(WC( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),  &  !*Vector
                          & -nGhosts,nElemsX+nGhosts+1, &  !*X bounds vector  
                          & -nGhosts,nElemsY+nGhosts+1,   &  !*Y bounds vector
                          & 1,nElemsX, &                   !*X bounds loop    
                          & 1,nElemsY)                       !*Y bounds loop 
#endif

#ifdef RECONSTRUCTFROMCONSERVED
CALL ReconstructionX_INPUT_CONSERVED(U( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&    !*Y bounds vector  
                          & 0,nElemsX+1, &                 !*X bounds loop  
                          & 1,nElemsY)                   !*Y bounds loop    

#ifdef PATHCONSERVATIVESHOCKDETECTION
CALL ReconstructionFixX_INPUT_CONSERVED(U( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&    !*Y bounds vector  
                          & 0,nElemsX+1, &                 !*X bounds loop  
                          & 1,nElemsY, &                   !*Y bounds loop    
                          & Cell_To_Limit_PCSD( -nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) )
#endif

CALL ConsToPrimInTheWholeMeshInterfaces(WM(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1),0,nElemsX+1+1,0,nElemsY+1+1,WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1))
WM(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)=WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)
CALL ConsToPrimInTheWholeMeshInterfaces(WP(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1),0,nElemsX+1+1,0,nElemsY+1+1,WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1))
WP(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)=WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)

#else
CALL ReconstructionX_INPUT_PRIMITIVE(WC( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&    !*Y bounds vector  
                          & 0,nElemsX+1, &                 !*X bounds loop  
                          & 1,nElemsY)                   !*Y bounds loop    

#ifdef PATHCONSERVATIVESHOCKDETECTION
CALL ReconstructionFixX_INPUT_PRIMITIVE(WC( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&    !*Y bounds vector  
                          & 0,nElemsX+1, &                 !*X bounds loop  
                          & 1,nElemsY, &                   !*Y bounds loop    
                          & Cell_To_Limit_PCSD( -nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) )
#endif

#endif

CALL NumericalFluxFX()

#ifdef RECONSTRUCTFROMCONSERVED
CALL ReconstructionY_INPUT_CONSERVED(U( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&  !*Y bounds vector 
                          & 1,nElemsX, &                   !*X bounds loop
                          & 0,nElemsY+1)                 !*Y bounds loop    


#ifdef PATHCONSERVATIVESHOCKDETECTION
CALL ReconstructionFixY_INPUT_CONSERVED(U( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&  !*Y bounds vector 
                          & 1,nElemsX, &                   !*X bounds loop
                          & 0,nElemsY+1, &                 !*Y bounds loop    
                          & Cell_To_Limit_PCSD( -nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) )
#endif


CALL ConsToPrimInTheWholeMeshInterfaces(WM(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1),0,nElemsX+1+1,0,nElemsY+1+1,WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1))
WM(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)=WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)
CALL ConsToPrimInTheWholeMeshInterfaces(WP(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1),0,nElemsX+1+1,0,nElemsY+1+1,WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1))
WP(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)=WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)


#else
CALL ReconstructionY_INPUT_PRIMITIVE(WC( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&  !*Y bounds vector 
                          & 1,nElemsX, &                   !*X bounds loop
                          & 0,nElemsY+1)                 !*Y bounds loop    

#ifdef PATHCONSERVATIVESHOCKDETECTION
CALL ReconstructionFixY_INPUT_PRIMITIVE(WC( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&  !*Y bounds vector 
                          & 1,nElemsX, &                   !*X bounds loop
                          & 0,nElemsY+1, &                 !*Y bounds loop    
                          & Cell_To_Limit_PCSD( -nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1) )
#endif

#endif



CALL NumericalFluxFY()

CALL SourceTerm_Conserved_INPUT_CONSERVED(t,U(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&    !*Y bounds vector  
                          & 1,nElemsX, &                   !*X bounds loop  
                          & 1,nElemsY, &                   !*Y bounds loop    
                          & MeshBary(1:nDims,1:nElemsX,1:nElemsY), &
                          & 1,nElemsX, &
                          & 1,nElemsY, &
                          & MeshGP   (1:nDims,1:nElemsX,1:nElemsY,1:nGPs,1:nGPs), &
                          & 1,nElemsX, &
                          & 1,nElemsY)

CALL UpdateTimeDerivative()


#elif defined(ACTIVEFLUX)


!*-------------------------
!*U
!*-------------------------
!*From cell averages W_X, W_Y to W_qp
!*FluxFX
!*FluxFY
!*SourceTerms
!*UpdateTimeDerivative(U)
!*-------------------------
#ifndef PRIMITIVEONLY

CALL FromCellAveragesToPointValues()
CALL FluxFX()                !*Filling FX for the update of U through reconstructed values
CALL FluxFY()                !*Filling FY for the update of U through reconstructed values
CALL SourceTerms(t)          !*Filling S for the update of U. Standard computation.
CALL UpdateTimeDerivative()  !*Filling Ut

#endif


#else
!*-------------------------
!*RMK: U(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
!*-------------------------
CALL ComputeTransitionMatricesConservedInput(U( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),  &  !*Vector
                          & -nGhosts,nElemsX+nGhosts+1, &  !*X bounds vector  
                          & -nGhosts,nElemsY+nGhosts+1,   &  !*Y bounds vector
                          & 1,nElemsX, &                   !*X bounds loop    
                          & 1,nElemsY)                       !*Y bounds loop 

#if(1==1)
CALL ReconstructionX_INPUT_CONSERVED(U( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&    !*Y bounds vector  
                          & 0,nElemsX+1, &                 !*X bounds loop  
                          & 1,nElemsY)                   !*Y bounds loop    
#else
CALL ShocksIndicatorX()
CALL ReconstructionX()
CALL ReconstructionFixX()
CALL PositivityLimiterX()
#endif
CALL ConsToPrimInTheWholeMeshInterfaces(WM(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1),0,nElemsX+1+1,0,nElemsY+1+1,WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1))
WM(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)=WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)
CALL ConsToPrimInTheWholeMeshInterfaces(WP(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1),0,nElemsX+1+1,0,nElemsY+1+1,WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1))
WP(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)=WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)


CALL NumericalFluxFX()

#if(1==1)
CALL ReconstructionY_INPUT_CONSERVED(U( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&  !*Y bounds vector 
                          & 1,nElemsX, &                   !*X bounds loop
                          & 0,nElemsY+1)                 !*Y bounds loop    
#else
CALL ShocksIndicatorY() 
CALL ReconstructionY()
CALL ReconstructionFixY()
CALL PositivityLimiterY()
#endif
CALL ConsToPrimInTheWholeMeshInterfaces(WM(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1),0,nElemsX+1+1,0,nElemsY+1+1,WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1))
WM(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)=WM_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)
CALL ConsToPrimInTheWholeMeshInterfaces(WP(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1),0,nElemsX+1+1,0,nElemsY+1+1,WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1))
WP(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)=WP_support(1:nVar,1:nGPs,0:nElemsX+1+1,0:nElemsY+1+1)


CALL NumericalFluxFY()

CALL SourceTerms(t)
CALL UpdateTimeDerivative()

#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE FVTimeDerivativeConservedOnly
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
SUBROUTINE FVTimeDerivativePrimitiveOnly(t,ReconstructionInput)
!-------------------------------------------------------------------------------!
#ifdef CENTEREDPRIMITIVE
USE MOD_Equation,       ONLY: BoundaryConditions_On_Cell_Centered_Data_Primitive_INPUT
USE MOD_FiniteVolume2D_vars, ONLY: WC
USE MOD_FiniteVolume2D_vars, ONLY: WCt
#endif

USE MOD_Equation,       ONLY: SourceTerm_Primitive_INPUT_PRIMITIVE

#ifdef ACTIVEFLUX
USE MOD_Equation,       ONLY: Impose_BC_on_W
USE MOD_FiniteVolume2D_vars, ONLY: W_X
USE MOD_FiniteVolume2D_vars, ONLY: W_Y
USE MOD_FiniteVolume2D_vars, ONLY: Wt_X
USE MOD_FiniteVolume2D_vars, ONLY: Wt_Y
#endif

!*Unified subroutines
USE MOD_Reconstruction,      ONLY: ReconstructionX_INPUT_PRIMITIVE
USE MOD_Reconstruction,      ONLY: ReconstructionY_INPUT_PRIMITIVE
USE MOD_FiniteVolume2D_vars, ONLY: nGhosts
USE MOD_FiniteVolume2D_vars, ONLY: nElemsX
USE MOD_FiniteVolume2D_vars, ONLY: nElemsY
USE MOD_FiniteVolume2D_vars, ONLY: U
USE MOD_FiniteVolume2D_vars, ONLY: nVar
USE MOD_FiniteVolume2D_vars, ONLY: Reconstruction

USE MOD_FiniteVolume2D_vars, ONLY: MeshBary
USE MOD_FiniteVolume2D_vars, ONLY: MeshGP
USE MOD_FiniteVolume2D_vars, ONLY: nDims
USE MOD_FiniteVolume2D_vars, ONLY: nGPs
#ifdef ACTIVEFLUX
USE MOD_FiniteVolume2D_vars, ONLY: MeshBary_X
USE MOD_FiniteVolume2D_vars, ONLY: MeshGP_X
USE MOD_FiniteVolume2D_vars, ONLY: MeshBary_Y
USE MOD_FiniteVolume2D_vars, ONLY: MeshGP_Y
#endif

!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,         INTENT(IN) :: t
INTEGER,      INTENT(IN),OPTIONAL :: ReconstructionInput
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ReconstructionUsed
#ifdef ACTIVEFLUX
REAL            :: W_X_support( 1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)!*NB:+1 in X direction with respect to U
REAL            :: W_Y_support( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)!*NB:+1 in Y direction with respect to U
#endif
INTEGER         :: ii,jj

#if defined(CENTEREDPRIMITIVE)

IF (PRESENT(ReconstructionInput)) THEN
  ReconstructionUsed=ReconstructionInput
ELSE
  ReconstructionUsed=Reconstruction
END IF


CALL BoundaryConditions_On_Cell_Centered_Data_Primitive_INPUT(WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),t)



!*=========================
!*WC
!*=========================
CALL ReconstructionX_INPUT_PRIMITIVE(WC( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&    !*Y bounds vector  
                          & 0,nElemsX+1, &                 !*X bounds loop  
                          & 1,nElemsY, &                   !*Y bounds loop    
                          & ReconstructionInput=ReconstructionUsed)

CALL NumericalFluxFX_W_INPUT( 0,nElemsX,&                    !*X bounds loop    
                            & 1,nElemsY)                     !*Y bounds loop 

CALL BSurfX_W_INPUT( 0,nElemsX,&                             !*X bounds loop 
                   & 1,nElemsY)                              !*Y bounds loop 

CALL BVolX_W_INPUT(WC( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                                    !*Vector
                  & -nGhosts,nElemsX+nGhosts+1, &            !*X bounds vector  
                  & -nGhosts,nElemsY+nGhosts+1,   &          !*Y bounds vector
                  & 1,nElemsX, &                             !*X bounds loop    
                  & 1,nElemsY, &                             !*Y bounds loop  
                  & ReconstructionInput=ReconstructionUsed)

CALL ReconstructionY_INPUT_PRIMITIVE(WC( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&  !*Y bounds vector 
                          & 1,nElemsX, &                   !*X bounds loop
                          & 0,nElemsY+1, &                 !*Y bounds loop    
                          & ReconstructionInput=ReconstructionUsed)


CALL NumericalFluxFY_W_INPUT( 1,nElemsX, &                   !*X bounds loop    
                            & 0,nElemsY )                    !*Y bounds loop  

CALL BSurfY_W_INPUT( 1,nElemsX,&                             !*X bounds loop 
                   & 0,nElemsY)                              !*Y bounds loop 

CALL BVolY_W_INPUT(WC( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                                    !*Vector
                  & -nGhosts,nElemsX+nGhosts+1, &            !*X bounds vector  
                  & -nGhosts,nElemsY+nGhosts+1,   &          !*Y bounds vector
                  & 1,nElemsX, &                             !*X bounds loop    
                  & 1,nElemsY, &                               !*Y bounds loop  
                  & ReconstructionInput=ReconstructionUsed)


CALL SourceTerm_Primitive_INPUT_PRIMITIVE(t,WC(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1,&    !*Y bounds vector  
                          & 1,nElemsX, &                   !*X bounds loop  
                          & 1,nElemsY, &                   !*Y bounds loop    
                          & MeshBary(1:nDims,1:nElemsX,1:nElemsY), &
                          & 1,nElemsX, &
                          & 1,nElemsY, &
                          & MeshGP   (1:nDims,1:nElemsX,1:nElemsY,1:nGPs,1:nGPs), &
                          & 1,nElemsX, &
                          & 1,nElemsY)

CALL UpdateTimeDerivative_W_INPUT(WCt(1:nVar,1:nElemsX,1:nElemsY), &    !*Vector
                                 & 1,nElemsX, &                         !*X bounds loop    
                                 & 1,nElemsY)                           !*Y bounds loop  






#elif defined(ACTIVEFLUX)
CALL Impose_BC_on_W(W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1),t)

W_X_support=W_X
W_Y_support=W_Y
#ifdef MOMENTUMINPRIMITIVEVARIABLES
#ifdef RECONSTRUCTIONVELOCITY
DO jj=-nGhosts,nElemsY+nGhosts+1
  DO ii=-nGhosts,nElemsX+nGhosts+1+1
    IF (W_X_support(1,ii,jj) .NE. 0.0) THEN
      W_X_support(2:3,ii,jj)=W_X_support(2:3,ii,jj)/W_X_support(1,ii,jj)
    END IF
  END DO
END DO
DO jj=-nGhosts,nElemsY+nGhosts+1+1
  DO ii=-nGhosts,nElemsX+nGhosts+1
    IF (W_Y_support(1,ii,jj) .NE. 0.0) THEN
      W_Y_support(2:3,ii,jj)=W_Y_support(2:3,ii,jj)/W_Y_support(1,ii,jj)
    END IF
  END DO
END DO
#endif
#endif

#ifdef MULTIFLUID
#ifdef RECONSTRUCTIONPHI
DO jj=-nGhosts,nElemsY+nGhosts+1
  DO ii=-nGhosts,nElemsX+nGhosts+1+1
    IF (W_X_support(1,ii,jj) .NE. 0.0) THEN
      W_X_support(nVar,ii,jj)=W_X_support(nVar,ii,jj)/W_X_support(1,ii,jj)
    END IF
  END DO
END DO
DO jj=-nGhosts,nElemsY+nGhosts+1+1
  DO ii=-nGhosts,nElemsX+nGhosts+1
    IF (W_Y_support(1,ii,jj) .NE. 0.0) THEN
      W_Y_support(nVar,ii,jj)=W_Y_support(nVar,ii,jj)/W_Y_support(1,ii,jj)
    END IF
  END DO
END DO
#endif
#endif



!*-------------------------
!*W_X
!*RMK: W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1)
!*-------------------------
!*ReconstructionX_W_X
!*NumericalFluxFX_W_X
!*BSurfX_W_X
!*BVolX_W_X
!*---
!*ReconstructionY_W_X
!*NumericalFluxFY_W_X
!*BSurfY_W_X
!*BVolY_W_X
!*---
!*UpdateTimeDerivative_W_X
!*-------------------------
CALL ComputeTransitionMatricesPrimitiveInput(W_X( 1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1),  &  !*Vector
                          & -nGhosts,nElemsX+nGhosts+1+1, &  !*X bounds vector  !*NB:+1 in X direction
                          & -nGhosts,nElemsY+nGhosts+1,   &  !*Y bounds vector
                          & 1,nElemsX+1, &                   !*X bounds loop    !*NB:+1 in X direction
                          & 1,nElemsY)                       !*Y bounds loop 



CALL ReconstructionX_INPUT_PRIMITIVE(W_X_support( 1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1), &                            !*Vector
                          & -nGhosts,nElemsX+nGhosts+1+1, &  !*X bounds vector  !*NB:+1 in X direction
                          & -nGhosts,nElemsY+nGhosts+1,   &  !*Y bounds vector
                          & 0,nElemsX+1+1, &                 !*X bounds loop    !*NB:+1 in X direction
                          & 1,nElemsY)                       !*Y bounds loop  

CALL Adjust_Interface_values()

CALL NumericalFluxFX_W_INPUT( 0,nElemsX+1,&                  !*X bounds loop    !*NB:+1 in X direction
                            & 1,nElemsY)                     !*Y bounds loop 


CALL BSurfX_W_INPUT( 0,nElemsX+1,&                           !*X bounds loop !*NB:+1 in X direction, like NumericalFluxFX_W_X
                   & 1,nElemsY)                              !*Y bounds loop 


CALL BVolX_W_INPUT(W_X( 1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1), &                                    !*Vector
                  & -nGhosts,nElemsX+nGhosts+1+1, &          !*X bounds vector  !*NB:+1 in X direction
                  & -nGhosts,nElemsY+nGhosts+1,   &          !*Y bounds vector
                  & 1,nElemsX+1, &                           !*X bounds loop    !*NB:+1 in X direction
                  & 1,nElemsY)                               !*Y bounds loop  

CALL ReconstructionY_INPUT_PRIMITIVE(W_X_support( 1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1), &                            !*Vector
                          & -nGhosts,nElemsX+nGhosts+1+1, &  !*X bounds vector  !*NB:+1 in X direction
                          & -nGhosts,nElemsY+nGhosts+1,   &  !*Y bounds vector
                          & 1,nElemsX+1, &                   !*X bounds loop    !*NB:+1 in X direction
                          & 0,nElemsY+1)                     !*Y bounds loop  

CALL Adjust_Interface_values()

CALL NumericalFluxFY_W_INPUT( 1,nElemsX+1, &                 !*X bounds loop    !*NB:+1 in X direction
                            & 0,nElemsY )                    !*Y bounds loop  

CALL BSurfY_W_INPUT( 1,nElemsX+1,&                           !*X bounds loop !*NB:+1 in X direction
                   & 0,nElemsY)                              !*Y bounds loop 

CALL BVolY_W_INPUT(W_X( 1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1), &                                    !*Vector
                  & -nGhosts,nElemsX+nGhosts+1+1, &          !*X bounds vector  !*NB:+1 in X direction
                  & -nGhosts,nElemsY+nGhosts+1,   &          !*Y bounds vector
                  & 1,nElemsX+1, &                           !*X bounds loop    !*NB:+1 in X direction
                  & 1,nElemsY)                               !*Y bounds loop  

CALL SourceTerm_Primitive_INPUT_PRIMITIVE(t,W_X(1:nVar,-nGhosts:nElemsX+nGhosts+1+1,-nGhosts:nElemsY+nGhosts+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1+1,  &  !*X bounds vector !*NB:+1 in X direction
                          &-nGhosts,nElemsY+nGhosts+1,&    !*Y bounds vector  
                          & 1,nElemsX+1, &                   !*X bounds loop !*NB:+1 in X direction
                          & 1,nElemsY, &                   !*Y bounds loop    
                          & MeshBary_X(1:nDims,1:nElemsX+1,1:nElemsY), &
                          & 1,nElemsX+1, &                  !*NB:+1 in X direction
                          & 1,nElemsY, &
                          & MeshGP_X(1:nDims,1:nElemsX+1,1:nElemsY,1:nGPs,1:nGPs), &
                          & 1,nElemsX+1, &                   !*NB:+1 in X direction
                          & 1,nElemsY)

!*Filling Wt_X
CALL UpdateTimeDerivative_W_INPUT(Wt_X(1:nVar,1:nElemsX+1,1:nElemsY), &   !*Vector
                                 & 1,nElemsX+1, &                         !*X bounds loop    !*NB:+1 in X direction
                                 & 1,nElemsY)                             !*Y bounds loop  
                                  

!*-------------------------
!*W_Y
!*RMK: W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1)
!*-------------------------
!*ReconstructionX_W_Y
!*NumericalFluxFX_W_Y
!*BSurfX_W_Y
!*BVolX_W_Y
!*---
!*ReconstructionY_W_Y
!*NumericalFluxFY_W_Y
!*BSurfY_W_Y
!*BVolY_W_Y
!*---
!*UpdateTimeDerivative_W_Y
!*-------------------------
CALL ComputeTransitionMatricesPrimitiveInput(W_Y( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1+1,&  !*Y bounds vector  !*NB:+1 in Y direction
                          & 1,nElemsX, &                 !*X bounds loop  
                          & 1,nElemsY+1)                 !*Y bounds loop    !*NB:+1 in Y direction


CALL ReconstructionX_INPUT_PRIMITIVE(W_Y_support( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1+1,&  !*Y bounds vector  !*NB:+1 in Y direction
                          & 0,nElemsX+1, &                 !*X bounds loop  
                          & 1,nElemsY+1)                   !*Y bounds loop    !*NB:+1 in Y direction

CALL Adjust_Interface_values()

CALL NumericalFluxFX_W_INPUT( 0,nElemsX,&                  !*X bounds loop    
                            & 1,nElemsY+1)                 !*Y bounds loop    !*NB:+1 in Y direction


CALL BSurfX_W_INPUT( 0,nElemsX,&                           !*X bounds loop 
                   & 1,nElemsY+1)                          !*Y bounds loop !*NB:+1 in Y direction, like NumericalFluxFX_W_Y

CALL BVolX_W_INPUT(W_Y( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1), &                                  !*Vector
                  &-nGhosts,nElemsX+nGhosts+1,  &          !*X bounds vector
                  &-nGhosts,nElemsY+nGhosts+1+1,&          !*Y bounds vector  !*NB:+1 in Y direction
                  & 1,nElemsX, &                           !*X bounds loop  
                  & 1,nElemsY+1)                           !*Y bounds loop    !*NB:+1 in Y direction


CALL ReconstructionY_INPUT_PRIMITIVE(W_Y_support( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1+1,&  !*Y bounds vector  !*NB:+1 in Y direction
                          & 1,nElemsX, &                   !*X bounds loop
                          & 0,nElemsY+1+1)                 !*Y bounds loop    !*NB:+1 in Y direction

CALL Adjust_Interface_values()

CALL NumericalFluxFY_W_INPUT( 1,nElemsX, &                 !*X bounds loop
                            & 0,nElemsY+1 )                !*Y bounds loop    !*NB:+1 in Y direction  


CALL BSurfY_W_INPUT( 1,nElemsX, &                          !*X bounds loop
                   & 0,nElemsY+1)                          !*Y bounds loop    !*NB:+1 in Y direction


CALL BVolY_W_INPUT(W_Y( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1), &                                  !*Vector
                  &-nGhosts,nElemsX+nGhosts+1,  &          !*X bounds vector
                  &-nGhosts,nElemsY+nGhosts+1+1,&          !*Y bounds vector  !*NB:+1 in Y direction
                  & 1,nElemsX, &                           !*X bounds loop  
                  & 1,nElemsY+1)                           !*Y bounds loop    !*NB:+1 in Y direction

CALL SourceTerm_Primitive_INPUT_PRIMITIVE(t,W_Y(1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1+1), &                          !*Vector
                          &-nGhosts,nElemsX+nGhosts+1,  &  !*X bounds vector
                          &-nGhosts,nElemsY+nGhosts+1+1,&    !*Y bounds vector !*NB:+1 in Y direction  
                          & 1,nElemsX, &                   !*X bounds loop
                          & 1,nElemsY+1, &                   !*Y bounds loop !*NB:+1 in Y direction    
                          & MeshBary_Y(1:nDims,1:nElemsX,1:nElemsY+1), &
                          & 1,nElemsX, &                  
                          & 1,nElemsY+1, &                  !*NB:+1 in Y direction
                          & MeshGP_Y(1:nDims,1:nElemsX,1:nElemsY+1,1:nGPs,1:nGPs), &
                          & 1,nElemsX, &                   
                          & 1,nElemsY+1)                      !*NB:+1 in Y direction


!*<========  !*Filling Wt_Y
CALL UpdateTimeDerivative_W_INPUT(Wt_Y(1:nVar,1:nElemsX,1:nElemsY+1), &  !*Vector
                                 & 1,nElemsX, &                          !*X bounds loop
                                 & 1,nElemsY+1)                          !*Y bounds loop    !*NB:+1 in Y direction  

#endif


!-------------------------------------------------------------------------------!
END SUBROUTINE FVTimeDerivativePrimitiveOnly
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
SUBROUTINE Identify_Interface_Cells()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_U
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_W_X
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_W_Y
USE MOD_FiniteVolume2D_vars,ONLY: Which_Fluid_In_Cell_U
USE MOD_FiniteVolume2D_vars,ONLY: N_Troubled_Neighbors
USE MOD_Equation           ,ONLY: Get_Gamma
USE MOD_Equation           ,ONLY: Impose_BC_on_Troubled_Cell_U
USE MOD_Equation           ,ONLY: Impose_BC_on_Troubled_Cell_W
USE MOD_Equation           ,ONLY: Impose_BC_on_Fluid_Cell_U
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj, indn
REAL               :: gmmL, gmmR, gmmD, gmmU
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

Troubled_Cell_U=0
Troubled_Cell_W_X=0
Troubled_Cell_W_Y=0

Which_Fluid_In_Cell_U=0

!*I add one layer
DO jj=0,nElemsY+1
  DO ii=0,nElemsX+1

      !*L
      gmmL=Get_Gamma(W_X(nVar,ii,jj))

      !*R
      gmmR=Get_Gamma(W_X(nVar,ii+1,jj))

      !*D
      gmmD=Get_Gamma(W_Y(nVar,ii,jj))

      !*U
      gmmU=Get_Gamma(W_Y(nVar,ii,jj+1))

      IF (  (gmmL .NE. gmmR) .OR. (gmmD .NE. gmmU) .OR. (gmmL .NE. gmmD) ) THEN !*<=== MIXED CELL
        Troubled_Cell_U(ii,jj)=1
        Which_Fluid_In_Cell_U(ii,jj)=2*INT( SIGN( 1.0, W_X(nVar,ii,jj)+W_X(nVar,ii+1,jj)+W_Y(nVar,ii,jj)+W_Y(nVar,ii,jj+1) ) )
        IF (Which_Fluid_In_Cell_U(ii,jj) .EQ. 0) THEN
          Which_Fluid_In_Cell_U(ii,jj)=2
        END IF
        !*====================================================
        ! IF ( ABS( Which_Fluid_In_Cell_U(ii,jj) ) .NE. 2 ) THEN
        !   PRINT*, "PROBLEM IN MIXED CELL"
        !   STOP
        ! END IF
        !*====================================================

      ELSE !*<=== SINGLE FLUID CELL

#ifndef LEVELSETINU
        !*SO THAT IT HAS THE CORRECT SIGN WHEN USED IN THE POST PROCESSING
        !*ONLY IF LEVELSET is not evolved in U
        U(nVar,ii,jj)=0.25*(W_X(nVar,ii,jj)+W_X(nVar,ii+1,jj)+W_Y(nVar,ii,jj)+W_Y(nVar,ii,jj+1))
#endif
        !*Otherwise it already has a sign

        Which_Fluid_In_Cell_U(ii,jj)=INT( SIGN(1.0,W_X(nVar,ii,jj)) ) !*I randomly take the sign by one of the surrounding primitive variables, it is the same
        IF (Which_Fluid_In_Cell_U(ii,jj) .EQ. 0) THEN
          Which_Fluid_In_Cell_U(ii,jj)=1
        END IF
        !*====================================================
        ! IF ( ABS( Which_Fluid_In_Cell_U(ii,jj) ) .NE. 1 ) THEN
        !   PRINT*, "PROBLEM IN SINGLE FLUID CELL"
        !   STOP
        ! END IF
        !*====================================================


      END IF

  END DO
END DO

!*================================
!*Detection of neighbours
!*================================
DO indn=1,N_Troubled_Neighbors
  CALL Impose_BC_on_Troubled_Cell_U()
  DO jj=0,nElemsY+1
    DO ii=0,nElemsX+1
      IF ( Troubled_Cell_U(ii,jj) .EQ. indn ) THEN !* Check neighbors
        IF ( Troubled_Cell_U(ii-1,jj) .EQ. 0 ) THEN
          Troubled_Cell_U(ii-1,jj)=indn+1
        END IF

        IF ( Troubled_Cell_U(ii+1,jj) .EQ. 0 ) THEN
          Troubled_Cell_U(ii+1,jj)=indn+1
        END IF

        IF ( Troubled_Cell_U(ii,jj-1) .EQ. 0 ) THEN
          Troubled_Cell_U(ii,jj-1)=indn+1
        END IF

        IF ( Troubled_Cell_U(ii,jj+1) .EQ. 0 ) THEN
          Troubled_Cell_U(ii,jj+1)=indn+1
        END IF

      END IF
    END DO
  END DO

END DO

CALL Impose_BC_on_Troubled_Cell_U()

!*================================
!*Detection of troubled in W_X and W_Y
!*================================

!*W_X
DO jj=1,nElemsY
  DO ii=1,nElemsX+1 !*NB:+1 in X direction

    IF ( (Troubled_Cell_U(ii-1,jj) .NE. 0) .OR. (Troubled_Cell_U(ii,jj) .NE. 0) ) THEN
      Troubled_Cell_W_X(ii,jj)=1
    END IF

  END DO
END DO


!*W_Y
DO jj=1,nElemsY+1 !*NB:+1 in Y direction
  DO ii=1,nElemsX

    IF ( (Troubled_Cell_U(ii,jj-1) .NE. 0) .OR. (Troubled_Cell_U(ii,jj) .NE. 0) ) THEN
      Troubled_Cell_W_Y(ii,jj)=1
    END IF

  END DO
END DO

CALL Impose_BC_on_Troubled_Cell_W()

CALL Impose_BC_on_Fluid_Cell_U()


! jj=2
! PRINT*
! DO ii=1,nElemsX
!   PRINT*, ii, Which_Fluid_In_Cell_U(ii,:)
! END DO


!-------------------------------------------------------------------------------!
END SUBROUTINE Identify_Interface_Cells
#endif
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
#ifdef MULTIFLUID
SUBROUTINE Replace_Conserved_At_Interface()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: Troubled_Cell_U
USE MOD_Equation,ONLY           : PrimToCons
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
REAL               :: Vtemp(nVar)
REAL               :: VtempL(nVar), VtempR(nVar), VtempD(nVar), VtempU(nVar)
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!


DO jj=1,nElemsY
  DO ii=1,nElemsX

      IF ( Troubled_Cell_U(ii,jj) .NE. 0 ) THEN
        CALL PrimToCons(W_X(1:nVar,ii,jj),VtempL)
        CALL PrimToCons(W_X(1:nVar,ii+1,jj),VtempR)
        CALL PrimToCons(W_Y(1:nVar,ii,jj),VtempD)
        CALL PrimToCons(W_Y(1:nVar,ii,jj+1),VtempU)
        U(1:5,ii,jj)=0.25*(VtempL(1:5)+VtempR(1:5)+VtempD(1:5)+VtempU(1:5))

        ! IF ( ABS(MeshBary(1,ii,jj)-0.6) .LT. 1e-2 ) THEN
        !   PRINT*, ii, jj
        !   PRINT*, VtempL
        !   PRINT*, VtempR
        !   PRINT*, VtempD
        !   PRINT*, VtempU
        !   PRINT*, U(1:5,ii,jj)
        ! END IF

      END IF

  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE Replace_Conserved_At_Interface
#endif
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef CENTEREDPRIMITIVE
SUBROUTINE Overwrite_WC_from_U()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructedVariable
USE MOD_FiniteVolume2D_vars,ONLY: LLX
USE MOD_FiniteVolume2D_vars,ONLY: RRX
USE MOD_FiniteVolume2D_vars,ONLY: LLY
USE MOD_FiniteVolume2D_vars,ONLY: RRY
USE MOD_FiniteVolume2D_vars,ONLY: MCtoP
USE MOD_FiniteVolume2D_vars,ONLY: MPtoC
USE MOD_FiniteVolume2D_vars,ONLY: t
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGP
USE MOD_Equation           ,ONLY: ConsToPrim
USE MOD_Equation           ,ONLY: BoundaryConditions
USE MOD_Reconstruction     ,ONLY: WENO_In_Quadrature_Points
#ifdef RELAXATION
USE MOD_FiniteVolume2D_vars,ONLY: EPS_LM
#endif
USE MOD_Equation           ,ONLY: Switching_Function
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj, indi, indj, iGP, jGP
REAL               :: U_qp(1:nVar,1:nGPs,1:nGPs) 
REAL               :: W_qp(1:nVar,1:nGPs,1:nGPs) 
REAL               :: WsupportBox(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts)
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: W_from_U(1:nVar)
REAL               :: theta
!-------------------------------------------------------------------------------!



CALL BoundaryConditions(t)

CALL ComputeTransitionMatricesConservedInput(U( 1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1),  &  !*Vector
                          & -nGhosts,nElemsX+nGhosts+1, &  !*X bounds vector  
                          & -nGhosts,nElemsY+nGhosts+1,   &  !*Y bounds vector
                          & 1,nElemsX, &                   !*X bounds loop    
                          & 1,nElemsY)                       !*Y bounds loop 

DO jj=1,nElemsY
  DO ii=1,nElemsX
    SELECT CASE(ABS(Reconstruction))
      CASE(1,10,2,20,21,22,23,24,25,26,27,28,29)
#if !defined(IMEX) || !defined(RELAXATION)
        CALL ConsToPrim(U(1:nVar,ii,jj),WC(1:nVar,ii,jj))
#else
        CALL ConsToPrim(U(1:nVar,ii,jj),W_from_U(1:nVar))
        theta=Switching_Function(EPS_LM)
        WC(1:nVar,ii,jj)=(1.0-theta)*W_from_U(1:nVar)+theta*WC(1:nVar,ii,jj)
#endif
      CASE(3,4,5,7)

#if(1==1)      
        PRINT*, "HO RECONSTRUCTION OF PRIMTIVE VARIABLES TO BE CHECKED"
        STOP
#endif

        SELECT CASE (ReconstructedVariable)
          CASE(0) !*Conserved
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=U(1:nVar,ii+indi,jj+indj)
              END DO
            END DO
          CASE(1) !*Characteristic
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=U(1:nVar,ii+indi,jj+indj)  !*Cons
              END DO
            END DO
          CASE(2) !*Primitive
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=MATMUL(MCtoP(:,:,ii,jj),U(1:nVar,ii+indi,jj+indj))  !*Cons->Prim
              END DO
            END DO
          CASE DEFAULT
            ErrorMessage = "Error in Overwrite_WC_from_U. ReconstructedVariable not available"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT

        CALL WENO_In_Quadrature_Points(&
                    WsupportBox(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts),& !*Entering with Cons
                    U_qp(1:nVar,1:nGPs,1:nGPs),&
                    Reconstruction,&
                    RRX(1:nVar,1:nVar,ii,jj),&
                    LLX(1:nVar,1:nVar,ii,jj),&
                    RRY(1:nVar,1:nVar,ii,jj),&
                    LLY(1:nVar,1:nVar,ii,jj))

        !*We transform back
        SELECT CASE (ReconstructedVariable)
          CASE(0) !*Conservative

            !*Do nothing

          CASE(1) !*Characteristic

            !*Do nothing

          CASE(2) !*Primitive

            DO jGP=1,nGPs
              DO iGP=1,nGPs
                U_qp(1:nVar,iGP,jGP)=MATMUL(MPtoC(:,:,ii,jj),U_qp(1:nVar,iGP,jGP))
              END DO
            END DO

          CASE DEFAULT
            ErrorMessage = "Error in Overwrite_WC_from_U. ReconstructedVariable not available"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT

        DO jGP=1,nGPS
          DO iGP=1,nGPS
            CALL ConsToPrim(U_qp(1:nVar,iGP,jGP),W_qp(1:nVar,iGP,jGP))
          END DO
        END DO
        DO jGP=1,nGPS
          DO iGP=1,nGPS
            WC(1:nVar,ii,jj)=WC(1:nVar,ii,jj)+WeightsGP(iGP,jGP)*W_qp(1:nVar,iGP,jGP)
          END DO
        END DO

      CASE DEFAULT
        ErrorMessage = "Error in Overwrite_WC_from_U. ReconstructedVariable not available"
        WRITE(*,*) ErrorMessage
        STOP

    END SELECT

  END DO !ii
END DO !jj

!-------------------------------------------------------------------------------!
END SUBROUTINE Overwrite_WC_from_U
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(CENTEREDPRIMITIVE) && defined(PATHCONSERVATIVESHOCKDETECTION)
SUBROUTINE Shock_Detector_Based_On_Path_Conservative()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: N_Limited_Neighbors_PCSD
USE MOD_FiniteVolume2D_vars,ONLY: K_coefficient_PCSD
USE MOD_FiniteVolume2D_vars,ONLY: Cell_To_Limit_PCSD
USE MOD_FiniteVolume2D_vars,ONLY: IsSomeoneFlagged_PCSD
USE MOD_FiniteVolume2D_vars,ONLY: Errors_PCSD
USE MOD_FiniteVolume2D_vars,ONLY: WhichVariableForLimiting_PCSD
USE MOD_Equation           ,ONLY: Impose_BC_on_Troubled_Cell_INPUT
USE MOD_Equation           ,ONLY: PrimToCons
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj, indn, indi, indj
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: WC_temp(1:nVar), Error_temp(1:nVar), Average_Error(1:nVar), weight
!-------------------------------------------------------------------------------!

IsSomeoneFlagged_PCSD=.FALSE.

Cell_To_Limit_PCSD=0
Average_Error(1:nVar)=0.0

CALL Compute_Errors_PCSD()

DO jj=1,nElemsY
  DO ii=1,nElemsX
    Average_Error(1:nVar)=Average_Error(1:nVar)+Errors_PCSD(1:nVar,ii,jj)
  END DO
END DO

Average_Error(1:nVar)=Average_Error(1:nVar)/(REAL(nElemsX)*REAL(nElemsY))


DO jj=1,nElemsY
  DO ii=1,nElemsX

    IF (WhichVariableForLimiting_PCSD .EQ. 1) THEN !*Density
      IF (Errors_PCSD(1,ii,jj) .GT. K_coefficient_PCSD*Average_Error(1)) THEN
        Cell_To_Limit_PCSD(ii,jj)=1
        IsSomeoneFlagged_PCSD=.TRUE.
      END IF
    ELSE IF (WhichVariableForLimiting_PCSD .EQ. 2) THEN !*Momentum
      IF ((Errors_PCSD(2,ii,jj) .GT. K_coefficient_PCSD*Average_Error(2)) .OR. (Errors_PCSD(3,ii,jj) .GT. K_coefficient_PCSD*Average_Error(3))) THEN !*NB: Checking momentum error
        Cell_To_Limit_PCSD(ii,jj)=1
        IsSomeoneFlagged_PCSD=.TRUE.
      END IF
    ELSE IF (WhichVariableForLimiting_PCSD .EQ. 3) THEN !*Energy
      IF (Errors_PCSD(4,ii,jj) .GT. K_coefficient_PCSD*Average_Error(4)) THEN
        Cell_To_Limit_PCSD(ii,jj)=1
        IsSomeoneFlagged_PCSD=.TRUE.
      END IF
    ELSE
      PRINT*, "Wrong choice for limiting variable."
      PRINT*, "Stopping in Shock_Detector_Based_On_Path_Conservative."
      PRINT*, WhichVariableForLimiting_PCSD
      STOP
    END IF

  END DO
END DO

!*================================
!*Detection of neighbours
!*================================
DO indn=1,N_Limited_Neighbors_PCSD
  CALL Impose_BC_on_Troubled_Cell_INPUT(Cell_To_Limit_PCSD(-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))
  DO jj=1,nElemsY
    DO ii=1,nElemsX
      IF ( Cell_To_Limit_PCSD(ii,jj) .EQ. indn ) THEN !* Check neighbors
        IF ( Cell_To_Limit_PCSD(ii-1,jj) .EQ. 0 ) THEN
          Cell_To_Limit_PCSD(ii-1,jj)=indn+1
        END IF

        IF ( Cell_To_Limit_PCSD(ii+1,jj) .EQ. 0 ) THEN
          Cell_To_Limit_PCSD(ii+1,jj)=indn+1
        END IF

        IF ( Cell_To_Limit_PCSD(ii,jj-1) .EQ. 0 ) THEN
          Cell_To_Limit_PCSD(ii,jj-1)=indn+1
        END IF

        IF ( Cell_To_Limit_PCSD(ii,jj+1) .EQ. 0 ) THEN
          Cell_To_Limit_PCSD(ii,jj+1)=indn+1
        END IF

      END IF
    END DO
  END DO

END DO

CALL Impose_BC_on_Troubled_Cell_INPUT(Cell_To_Limit_PCSD(-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1))



!-------------------------------------------------------------------------------!
END SUBROUTINE Shock_Detector_Based_On_Path_Conservative
#endif
!===============================================================================!
!
!
!
!===============================================================================!
#if defined(CENTEREDPRIMITIVE)
SUBROUTINE Compute_Errors_PCSD()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: WC
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: Errors_PCSD
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
USE MOD_Equation           ,ONLY: PrimToCons
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj, indn, indi, indj
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: WC_temp(1:nVar), Error_temp(1:nVar), weight
!-------------------------------------------------------------------------------!


Errors_PCSD=0.0

DO jj=1,nElemsY
  DO ii=1,nElemsX

    !*Error in Conserved variables
    Error_temp(1:nVar)=0.0

#ifndef SMOOTHINGLIMITER
    CALL PrimToCons(WC(1:nVar,ii,jj),WC_temp(1:nVar))
    Error_temp(1:nVar)=ABS(U(1:nVar,ii,jj)-WC_temp(1:nVar))
#else
    !*If we are at the boundary, only consider the cell itself, as corners may not be available
    IF ( (ii .EQ. 1) .OR. (ii .EQ. nElemsX) .OR. (jj .EQ. 1) .OR. (jj .EQ. nElemsY)  ) THEN
      CALL PrimToCons(WC(1:nVar,ii,jj),WC_temp(1:nVar))
      Error_temp(1:nVar)=ABS(U(1:nVar,ii,jj)-WC_temp(1:nVar))
    ELSE
      !*Cell itself
      indi=ii
      indj=jj
      weight=4.0/9.0
      CALL PrimToCons(WC(1:nVar,indi,indj),WC_temp(1:nVar))
      Error_temp(1:nVar)=Error_temp(1:nVar)+ABS(U(1:nVar,indi,indj)-WC_temp(1:nVar))*weight

      !*Left
      indi=ii-1
      indj=jj
      weight=1.0/9.0
      CALL PrimToCons(WC(1:nVar,indi,indj),WC_temp(1:nVar))
      Error_temp(1:nVar)=Error_temp(1:nVar)+ABS(U(1:nVar,indi,indj)-WC_temp(1:nVar))*weight

      !*Right
      indi=ii+1
      indj=jj
      weight=1.0/9.0
      CALL PrimToCons(WC(1:nVar,indi,indj),WC_temp(1:nVar))
      Error_temp(1:nVar)=Error_temp(1:nVar)+ABS(U(1:nVar,indi,indj)-WC_temp(1:nVar))*weight

      !*Bottom
      indi=ii
      indj=jj-1
      weight=1.0/9.0
      CALL PrimToCons(WC(1:nVar,indi,indj),WC_temp(1:nVar))
      Error_temp(1:nVar)=Error_temp(1:nVar)+ABS(U(1:nVar,indi,indj)-WC_temp(1:nVar))*weight

      !*Up
      indi=ii
      indj=jj+1
      weight=1.0/9.0
      CALL PrimToCons(WC(1:nVar,indi,indj),WC_temp(1:nVar))
      Error_temp(1:nVar)=Error_temp(1:nVar)+ABS(U(1:nVar,indi,indj)-WC_temp(1:nVar))*weight

      !*Bottom left
      indi=ii-1
      indj=jj-1
      weight=1.0/36.0
      CALL PrimToCons(WC(1:nVar,indi,indj),WC_temp(1:nVar))
      Error_temp(1:nVar)=Error_temp(1:nVar)+ABS(U(1:nVar,indi,indj)-WC_temp(1:nVar))*weight

      !*Bottom right
      indi=ii+1
      indj=jj-1
      weight=1.0/36.0
      CALL PrimToCons(WC(1:nVar,indi,indj),WC_temp(1:nVar))
      Error_temp(1:nVar)=Error_temp(1:nVar)+ABS(U(1:nVar,indi,indj)-WC_temp(1:nVar))*weight

      !*Up left
      indi=ii-1
      indj=jj+1
      weight=1.0/36.0
      CALL PrimToCons(WC(1:nVar,indi,indj),WC_temp(1:nVar))
      Error_temp(1:nVar)=Error_temp(1:nVar)+ABS(U(1:nVar,indi,indj)-WC_temp(1:nVar))*weight

      !*Up right
      indi=ii+1
      indj=jj+1
      weight=1.0/36.0
      CALL PrimToCons(WC(1:nVar,indi,indj),WC_temp(1:nVar))
      Error_temp(1:nVar)=Error_temp(1:nVar)+ABS(U(1:nVar,indi,indj)-WC_temp(1:nVar))*weight

    END IF
#endif

    Errors_PCSD(1:nVar,ii,jj)=Error_temp(1:nVar)

  END DO
END DO

#ifdef SMOOTHINGLIMITER
!*MAKE SURE A 1D TEST IS REALLY A 1D TEST
#if(1==1)
IF (   (InitialCondition .EQ. (-302)) .OR. (InitialCondition .EQ. (-303)) .OR. (InitialCondition .EQ. (350)) &
& .OR. (InitialCondition .EQ. (-304)) .OR. (InitialCondition .EQ. (-305)) ) THEN

#if(1==0)
PRINT*, "I prefer to keep it deactivated to be consistent with 2D"
PRINT*, "However, it did not solve the problem with left boundary of tests -302 and -304."
PRINT*, "Stop in finitevolume in Compute_Errors_PCSD"
STOP

DO jj=1,nElemsY

  Error_temp(1:nVar)=0.0

  !*===============
  !*LEFT BOUNDARY
  !*===============
  !*Cell itself
  indi=1
  indj=jj
  weight=2.0/3.0
  CALL PrimToCons(WC(1:nVar,indi,indj),WC_temp(1:nVar))
  Error_temp(1:nVar)=Error_temp(1:nVar)+ABS(U(1:nVar,indi,indj)-WC_temp(1:nVar))*weight

  !*No left

  !*Right
  indi=2
  indj=jj
  weight=1.0/6.0
  CALL PrimToCons(WC(1:nVar,indi,indj),WC_temp(1:nVar))
  Error_temp(1:nVar)=Error_temp(1:nVar)+ABS(U(1:nVar,indi,indj)-WC_temp(1:nVar))*weight

  Errors_PCSD(1:nVar,1,jj)=Error_temp(1:nVar)

  !*===============
  !*RIGHT BOUNDARY
  !*===============
  Error_temp(1:nVar)=0.0

  !*Cell itself
  indi=nElemsX
  indj=jj
  weight=2.0/3.0
  CALL PrimToCons(WC(1:nVar,indi,indj),WC_temp(1:nVar))
  Error_temp(1:nVar)=Error_temp(1:nVar)+ABS(U(1:nVar,indi,indj)-WC_temp(1:nVar))*weight

  !*No right

  !*Left
  indi=nElemsX-1
  indj=jj
  weight=1.0/6.0
  CALL PrimToCons(WC(1:nVar,indi,indj),WC_temp(1:nVar))
  Error_temp(1:nVar)=Error_temp(1:nVar)+ABS(U(1:nVar,indi,indj)-WC_temp(1:nVar))*weight

  Errors_PCSD(1:nVar,nElemsX,jj)=Error_temp(1:nVar)

END DO
#endif

  Errors_PCSD(1:nVar,-nGhosts:nElemsX+nGhosts+1,1)       = Errors_PCSD(1:nVar,-nGhosts:nElemsX+nGhosts+1,2)
  Errors_PCSD(1:nVar,-nGhosts:nElemsX+nGhosts+1,nElemsY) = Errors_PCSD(1:nVar,-nGhosts:nElemsX+nGhosts+1,nElemsY-1)

END IF
#endif
#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE Compute_Errors_PCSD
#endif
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE Impose_Symmetric_Update_Y_Axis()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
#endif
#ifdef ACTIVEFLUX
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
#endif
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: t
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
USE MOD_Equation           ,ONLY: BoundaryConditions
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj, indn, indi, indj
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: W_temp(1:nVar)
!-------------------------------------------------------------------------------!



W_temp(1:nVar)=0.0

DO jj=1,nElemsY
  DO ii=1,nElemsX/2

#ifdef CENTEREDPRIMITIVE

#endif

  END DO
END DO


#ifdef ACTIVEFLUX
!*W_X
DO jj=1,nElemsY
  DO ii=1,(nElemsX+1)/2


  END DO
END DO

!*W_Y
DO jj=1,nElemsY+1
  DO ii=1,nElemsX/2


  END DO
END DO
#endif

CALL BoundaryConditions(t)

!-------------------------------------------------------------------------------!
END SUBROUTINE Impose_Symmetric_Update_Y_Axis
!===============================================================================!
!
!
!
!===============================================================================!
END MODULE MOD_FiniteVolume2D
!-------------------------------------------------------------------------------!
