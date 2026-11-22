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



#ifdef CENTEREDPRIMITIVE
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
#ifdef CENTEREDPRIMITIVE
PUBLIC :: FVTimeDerivativePrimitiveOnly
PUBLIC :: FVTimeDerivativeConservedOnly
#endif
#if defined(CENTEREDPRIMITIVE)
PUBLIC :: Compute_Errors_PCSD
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

#ifdef CENTEREDPRIMITIVE
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





#ifdef CENTEREDPRIMITIVE
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




#ifdef CENTEREDPRIMITIVE
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
#ifdef CENTEREDPRIMITIVE
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


#ifdef CENTEREDPRIMITIVE
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


#ifdef CENTEREDPRIMITIVE
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
#ifdef CENTEREDPRIMITIVE
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
#ifdef CENTEREDPRIMITIVE
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
#ifdef CENTEREDPRIMITIVE
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
#ifdef CENTEREDPRIMITIVE
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
#ifdef CENTEREDPRIMITIVE
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
#ifdef CENTEREDPRIMITIVE
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
#ifdef CENTEREDPRIMITIVE
SUBROUTINE FVTimeDerivativeConservedOnly(t,Optional_Input)
!-------------------------------------------------------------------------------!
USE MOD_Equation,       ONLY: BoundaryConditions
USE MOD_Equation,       ONLY: SourceTerms
USE MOD_Equation,       ONLY: SourceTerm_Conserved_INPUT_CONSERVED
#ifdef CENTEREDPRIMITIVE
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
#ifdef CENTEREDPRIMITIVE
SUBROUTINE FVTimeDerivativePrimitiveOnly(t,ReconstructionInput)
!-------------------------------------------------------------------------------!
#ifdef CENTEREDPRIMITIVE
USE MOD_Equation,       ONLY: BoundaryConditions_On_Cell_Centered_Data_Primitive_INPUT
USE MOD_FiniteVolume2D_vars, ONLY: WC
USE MOD_FiniteVolume2D_vars, ONLY: WCt
#endif

USE MOD_Equation,       ONLY: SourceTerm_Primitive_INPUT_PRIMITIVE


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






#endif


!-------------------------------------------------------------------------------!
END SUBROUTINE FVTimeDerivativePrimitiveOnly
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
