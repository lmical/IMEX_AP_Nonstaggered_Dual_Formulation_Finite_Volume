!===============================================================================!
MODULE MOD_Mesh
!===============================================================================!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE BuildMesh
  MODULE PROCEDURE BuildMesh
END INTERFACE
INTERFACE GlobalElem
  MODULE PROCEDURE GlobalElem
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: BuildMesh
PUBLIC :: GlobalElem
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
SUBROUTINE BuildMesh()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: MESH_SX
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: NormVectX
USE MOD_FiniteVolume2D_vars,ONLY: NormVectY
USE MOD_FiniteVolume2D_vars,ONLY: TangVectX
USE MOD_FiniteVolume2D_vars,ONLY: TangVectY
USE MOD_FiniteVolume2D_vars,ONLY: MeshNodes
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP   
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGP   
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGPBnd  

!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! Local Variables
INTEGER :: ii, jj, iGP, jGP
REAL, DIMENSION(nGPs) :: quadWeights1D, quadNodes1D 
!-------------------------------------------------------------------------------!

MeshNodes = 0.0
MeshBary  = 0.0

Mesh_SX    = ABS(Mesh_X1-Mesh_X0)
Mesh_DX(1) = ABS(Mesh_SX(1))/(REAL(nElemsX))
Mesh_DX(2) = ABS(Mesh_SX(2))/(REAL(nElemsY))

DO jj=0,nElemsY
  DO ii=0,nElemsX
    MeshNodes(1:nDims,ii,jj) = Mesh_X0(1:2) + (/REAL(ii),REAL(jj)/)*Mesh_DX(1:2)
  END DO
END DO

DO jj=1,nElemsY
  DO ii=1,nElemsX
    MeshBary(1:nDims,ii,jj) = MeshNodes(1:nDims,ii-1,jj-1) + 0.5*Mesh_DX(1:2)
  END DO
END DO


!*------------------------------
!*NB: +1 in X and Y direction added for Active flux in all the normals
!*------------------------------
!------------------------------!
! Normal vectors: x-direction  !
!------------------------------!
NormVectX(1:nDims) = (/1.0,0.0/)


!------------------------------!
! Normal vector: y-direction  !
!------------------------------!
NormVectY(1:nDims) = (/0.0,1.0/)

!------------------------------!
! Tangent vectors: x-direction !
!------------------------------!
TangVectX(1:nDims) = (/0.0,1.0/)

!------------------------------!
! Tangent vectors: y-direction !
!------------------------------!
TangVectY(1:nDims) = (/-1.0,0.0/)

  !------------------------------!
  !   MeshGP Quadrature          !
  !------------------------------!

SELECT CASE (nGPs)
  CASE(1)
    quadWeights1D(1) = 1.0
    quadNodes1D(1)  =  0.0 
  CASE(2)
    quadWeights1D = (/0.5,0.5 /)
    quadNodes1D  =  (/- 1./(2.*sqrt(3.)), 1./(2.*sqrt(3.)) /)
  CASE(3) !*NOT OK FOR WENO5 because nonlinear weights arise
    quadWeights1D = (/5./18.,4./9., 5./18. /)
    quadNodes1D  =  (/- 0.5*sqrt(3./5.), 0.,  0.5*sqrt(3./5.) /)
  CASE(4)
    quadWeights1D = (/  (18.-SQRT(30.))/72., (18.+SQRT(30.))/72., (18.+SQRT(30.))/72. , (18.-SQRT(30.))/72. /)
    quadNodes1D  =  (/  -0.5*SQRT(3./7.+2./7.*SQRT(6./5.)), -0.5*SQRT(3./7.-2./7.*SQRT(6./5.)),&
     0.5*SQRT(3./7.-2./7.*SQRT(6./5.)), 0.5*SQRT(3./7.+2./7.*SQRT(6./5.))  /)
  CASE DEFAULT
    PRINT*, "Quadrature not implemented"
    STOP
END SELECT


DO iGP = 1, nGPs
  WeightsGPBnd(iGP) = quadWeights1D(iGP)
  DO jGP = 1, nGPs
    WeightsGP(iGP,jGP) = quadWeights1D(iGP)* quadWeights1D(jGP)
    DO jj=1,nElemsY
      DO ii=1,nElemsX
        MeshGP(1:nDims,ii,jj,iGP,jGP) = (/ MeshBary(1,ii,jj) +quadNodes1D(iGP)*Mesh_DX(1) , MeshBary(2,ii,jj)  +quadNodes1D(jGP)*Mesh_DX(2) /)
      END DO
    END DO
  END DO
END DO



!-------------------------------------------------------------------------------!
END SUBROUTINE BuildMesh
!===============================================================================!
!
!
!
!===============================================================================!
INTEGER FUNCTION GlobalElem(ii,jj)
!-------------------------------------------------------------------------------!

USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! Local Variables
INTEGER :: ii, jj
!-------------------------------------------------------------------------------!
GlobalElem = MODULO(jj-1,nElemsY)*NelemsX + MODULO(ii-1,NelemsX) + 1
!-------------------------------------------------------------------------------!
END FUNCTION GlobalElem
!===============================================================================!
!
!===============================================================================!
END MODULE MOD_Mesh
!-------------------------------------------------------------------------------!
