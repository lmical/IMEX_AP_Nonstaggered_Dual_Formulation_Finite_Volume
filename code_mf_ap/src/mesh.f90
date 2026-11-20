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
#ifdef ACTIVEFLUX
!*For staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary_X
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_X

!*For staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary_Y
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP_Y
#endif

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

#ifdef ACTIVEFLUX
!*For staggering in X direction
DO jj=1,nElemsY
  DO ii=1,nElemsX
    MeshBary_X(2,ii,jj) = MeshBary(2,ii,jj)                !*Same y
    MeshBary_X(1,ii,jj) = MeshBary(1,ii,jj)-0.5*Mesh_DX(1) !*Staggered x
  END DO
END DO

!*Last point in X direction skipped
DO jj=1,nElemsY
  MeshBary_X(2,nElemsX+1,jj) = MeshBary(2,nElemsX,jj)                !*Same y
  MeshBary_X(1,nElemsX+1,jj) = MeshBary(1,nElemsX,jj)+0.5*Mesh_DX(1) !*Staggered x
END DO

!*------------------------------
!*SAFETY CHECK
!*------------------------------
!*Same y
!*Staggered x
!*------------------------------
! DO jj=1,nElemsY
!   PRINT*
!   PRINT*, "jj", jj, MeshBary_X(2,1,jj) !*y level
!   DO ii=1,nElemsX+1
!     PRINT*, ii, jj, MeshBary_X(:,ii,jj)
!   END DO
! END DO
! STOP
!*------------------------------


!*For staggering in Y direction
DO jj=1,nElemsY
  DO ii=1,nElemsX
    MeshBary_Y(1,ii,jj) = MeshBary(1,ii,jj)                !*Same x
    MeshBary_Y(2,ii,jj) = MeshBary(2,ii,jj)-0.5*Mesh_DX(2) !*Staggered y
  END DO
END DO

!*Last point in Y direction skipped
DO ii=1,nElemsX
  MeshBary_Y(1,ii,nElemsY+1) = MeshBary(1,ii,nElemsY)                !*Same x
  MeshBary_Y(2,ii,nElemsY+1) = MeshBary(2,ii,nElemsY)+0.5*Mesh_DX(2) !*Staggered y
END DO

!*------------------------------
!*SAFETY CHECK
!*------------------------------
!*Same x
!*Staggered y
!*------------------------------
! DO jj=1,nElemsY+1
!   PRINT*
!   PRINT*, "jj", jj, MeshBary_Y(2,1,jj)
!   DO ii=1,nElemsX
!     PRINT*, ii, jj, MeshBary_Y(:,ii,jj)
!   END DO
! END DO
! STOP
!*------------------------------
#endif

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


#ifdef ACTIVEFLUX
!*For staggering in X direction
DO iGP = 1, nGPs
  DO jGP = 1, nGPs
    DO jj=1,nElemsY
      DO ii=1,nElemsX+1
        MeshGP_X(1:nDims,ii,jj,iGP,jGP) = (/ MeshBary_X(1,ii,jj) +quadNodes1D(iGP)*Mesh_DX(1) , MeshBary_X(2,ii,jj)  +quadNodes1D(jGP)*Mesh_DX(2) /)
        !*------------------------
        !*SAFETY CHECK
        !*------------------------
        ! IF (ii .NE. nElemsX+1) THEN
        !   IF( ABS( MeshGP(1,ii,jj,iGP,jGP) - MeshGP_X(1,ii,jj,iGP,jGP)-0.5*Mesh_DX(1) ) .GT. 1e-15 ) THEN
        !     PRINT*, "PROBLEM 1",ii,jj,iGP,jGP
        !     PRINT*, ABS(MeshGP(1,ii,jj,iGP,jGP) - MeshGP_X(1,ii,jj,iGP,jGP)-0.5*Mesh_DX(1))
        !     STOP
        !   END IF
        !   IF( ABS( MeshGP(2,ii,jj,iGP,jGP) - MeshGP_X(2,ii,jj,iGP,jGP) ) .GT. 1e-15 ) THEN
        !     PRINT*, "PROBLEM 2",ii,jj,iGP,jGP
        !     PRINT*, ABS( MeshGP(2,ii,jj,iGP,jGP) - MeshGP_X(2,ii,jj,iGP,jGP) )
        !     STOP
        !   END IF
        ! ELSE
        !   IF( ABS( MeshGP(1,nElemsX,jj,iGP,jGP) - MeshGP_X(1,nElemsX+1,jj,iGP,jGP)+0.5*Mesh_DX(1) ) .GT. 1e-15 ) THEN
        !     PRINT*, "PROBLEM 3",ii,jj,iGP,jGP
        !     PRINT*, ABS( MeshGP(1,nElemsX,jj,iGP,jGP) - MeshGP_X(1,nElemsX+1,jj,iGP,jGP)+0.5*Mesh_DX(1) )
        !     STOP
        !   END IF
        !   IF( ABS( MeshGP(2,nElemsX,jj,iGP,jGP) - MeshGP_X(2,nElemsX+1,jj,iGP,jGP) ) .GT. 1e-15 ) THEN
        !     PRINT*, "PROBLEM 4",ii,jj,iGP,jGP
        !     PRINT*, ABS( MeshGP(2,nElemsX,jj,iGP,jGP) - MeshGP_X(2,nElemsX+1,jj,iGP,jGP) )
        !     STOP
        !   END IF
        ! END IF
        !*------------------------
      END DO
    END DO
  END DO
END DO


!*For staggering in Y direction
DO iGP = 1, nGPs
  DO jGP = 1, nGPs
    DO jj=1,nElemsY+1
      DO ii=1,nElemsX
        MeshGP_Y(1:nDims,ii,jj,iGP,jGP) = (/ MeshBary_Y(1,ii,jj) +quadNodes1D(iGP)*Mesh_DX(1) , MeshBary_Y(2,ii,jj)  +quadNodes1D(jGP)*Mesh_DX(2) /)
        !*------------------------
        !*SAFETY CHECK
        !*------------------------
        ! IF (jj .NE. nElemsY+1) THEN
        !   IF( ABS( MeshGP(2,ii,jj,iGP,jGP) - MeshGP_Y(2,ii,jj,iGP,jGP)-0.5*Mesh_DX(2) ) .GT. 1e-15 ) THEN
        !     PRINT*, "PROBLEM 5",ii,jj,iGP,jGP
        !     PRINT*, ABS(MeshGP(2,ii,jj,iGP,jGP) - MeshGP_X(2,ii,jj,iGP,jGP)-0.5*Mesh_DX(2))
        !     STOP
        !   END IF
        !   IF( ABS( MeshGP(1,ii,jj,iGP,jGP) - MeshGP_Y(1,ii,jj,iGP,jGP) ) .GT. 1e-15 ) THEN
        !     PRINT*, "PROBLEM 6",ii,jj,iGP,jGP
        !     PRINT*, ABS( MeshGP(1,ii,jj,iGP,jGP) - MeshGP_Y(1,ii,jj,iGP,jGP) )
        !     STOP
        !   END IF
        ! ELSE
        !   IF( ABS( MeshGP(2,ii,nElemsY,iGP,jGP) - MeshGP_Y(2,ii,nElemsY+1,iGP,jGP)+0.5*Mesh_DX(2) ) .GT. 1e-15 ) THEN
        !     PRINT*, "PROBLEM 7",ii,jj,iGP,jGP
        !     PRINT*, ABS( MeshGP(2,ii,nElemsY,iGP,jGP) - MeshGP_Y(2,ii,nElemsY+1,iGP,jGP)+0.5*Mesh_DX(2) )
        !     STOP
        !   END IF
        !   IF( ABS( MeshGP(1,ii,nElemsY,iGP,jGP) - MeshGP_Y(1,ii,nElemsY+1,iGP,jGP) ) .GT. 1e-15 ) THEN
        !     PRINT*, "PROBLEM 8",ii,jj,iGP,jGP
        !     PRINT*, ABS( MeshGP(1,ii,nElemsY,iGP,jGP) - MeshGP_Y(1,ii,nElemsY+1,iGP,jGP) )
        !     STOP
        !   END IF
        ! END IF
        !*------------------------
      END DO
    END DO
  END DO
END DO



#endif

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
