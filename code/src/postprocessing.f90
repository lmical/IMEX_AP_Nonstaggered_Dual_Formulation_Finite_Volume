!===============================================================================!
MODULE MOD_PostProcessing
!===============================================================================!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE PostProcessing  
  MODULE PROCEDURE PostProcessing 
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: PostProcessing
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
SUBROUTINE PostProcessing
USE MOD_Equation, ONLY :  ExactFunction
USE MOD_Equation, ONLY :  ConsToPrim
#if defined(CENTEREDPRIMITIVE) && defined(PATHCONSERVATIVESHOCKDETECTION)
USE MOD_FiniteVolume2D, ONLY : Compute_Errors_PCSD
#endif
USE MOD_FiniteVolume2D_vars
IMPLICIT NONE
INTEGER              :: ii, jj, iGP, jGP, iVar, i, j, k
REAL, ALLOCATABLE    :: Uexact(:,:,:), L1error(:)
REAL                 :: Error_Max_temp(1:nVar)
CHARACTER(LEN=255)   :: FNAME
REAL, DIMENSION(nVar, nGPs, nGPs) :: Utemp
REAL                 :: averageJacobi, std_deviation

#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
REAL, ALLOCATABLE    :: Wexact(:,:,:)
REAL, DIMENSION(nVar, nGPs, nGPs) :: WtempCons, WtempPrim
#endif

PRINT*, "--------------------------"
#ifndef PRIMITIVEONLY
ALLOCATE( L1error(1:nVar) )
ALLOCATE( Uexact(1:nVar,1:nElemsX,1:nElemsY))
 
Uexact(:,:,:) = 0.
L1error(:) = 0. 
Utemp = 0.
DO jj = 1,nElemsY  
  DO ii = 1,nElemsX  
    DO iGP=1,nGPs
      DO jGP=1,nGPs
        CALL ExactFunction(InitialCondition, tEnd, MeshGP(:,ii,jj,iGP,jGP), Utemp(1:nVar,iGP,jGP))
        Uexact(1:nVar,ii,jj) = Uexact(1:nVar,ii,jj) + WeightsGP(iGP,jGP)* Utemp(1:nVar,iGP,jGP)
      END DO
    END DO
    L1error(:) = L1error(:) + ABS( U(:,ii,jj) - Uexact(:,ii,jj) ) 
  ENDDO
ENDDO

L1error(:) = L1error(:) * MESH_DX(1) * MESH_DX(2)

PRINT*, "Error", L1error


FNAME = 'ErrorL1_XXXX_XXXX.dat'
WRITE(FNAME(9:12),FMT="(I4.4)") nElemsX
WRITE(FNAME(14:17),FMT="(I4.4)") nElemsY

OPEN(666,FILE=TRIM(FNAME))
WRITE(666,*) nElemsX, nElemsY, L1error(1:nVar)

DEALLOCATE( L1error )
DEALLOCATE( Uexact )
#endif

!*One can uncomment them if needed. 
!*They work
! WRITE(666,*) 'computational time'
! WRITE(666,*) computationalTime
! WRITE(666,*) 'minimum ro'
! WRITE(666,*) global_min


#if defined(CENTEREDPRIMITIVE) && defined(PATHCONSERVATIVESHOCKDETECTION)
Error_Max_temp = 0.0

CALL Compute_Errors_PCSD()

DO iVar=1,nVar
  Error_Max_temp(iVar)=MAXVAL(Errors_PCSD(iVar,:,:))
END DO

FNAME = 'WC_MaxError_XXXX_XXXX.dat'
WRITE(FNAME(13:16),FMT="(I4.4)") nElemsX
WRITE(FNAME(18:21),FMT="(I4.4)") nElemsY

OPEN(666,FILE=TRIM(FNAME))
WRITE(666,*) nElemsX, nElemsY, Error_Max_temp(1:nVar)
#endif


#ifdef ACTIVEFLUX

!*Staggering in X direction
ALLOCATE( L1error(1:nVar) )
ALLOCATE( Wexact(1:nVar,1:nElemsX+1,1:nElemsY)) !*NB:+1 in X direction
 
Wexact(:,:,:) = 0.
L1error(:) = 0. 
WtempCons = 0.
WtempPrim = 0.
DO jj = 1,nElemsY  
  DO ii = 1,nElemsX+1 !*NB:+1 in X direction
    DO iGP=1,nGPs
      DO jGP=1,nGPs
        CALL ExactFunction(InitialCondition, tEnd, MeshGP_X(:,ii,jj,iGP,jGP), WtempCons(1:nVar,iGP,jGP))
        CALL ConsToPrim(WtempCons(1:nVar,iGP,jGP),WtempPrim(1:nVar,iGP,jGP))
        Wexact(1:nVar,ii,jj) = Wexact(1:nVar,ii,jj) + WeightsGP(iGP,jGP)* WtempPrim(1:nVar,iGP,jGP)
      END DO
    END DO
    L1error(:) = L1error(:) + ABS( W_X(:,ii,jj) - Wexact(:,ii,jj) ) 
  ENDDO
ENDDO

L1error(:) = L1error(:) * MESH_DX(1) * MESH_DX(2)

PRINT*, "Error W_X", L1error

FNAME = 'W_X_ErrorL1_XXXX_XXXX.dat'
WRITE(FNAME(9+4:12+4),FMT="(I4.4)") nElemsX
WRITE(FNAME(14+4:17+4),FMT="(I4.4)") nElemsY

OPEN(666,FILE=TRIM(FNAME))
WRITE(666,*) nElemsX, nElemsY, L1error(1:nVar)

#ifdef ACTIVEFLUX
#ifdef IMEX
WRITE(666,*) HyperbolicityLoss
#endif
#endif

DEALLOCATE( L1error )
DEALLOCATE( Wexact )


!*Staggering in Y direction
ALLOCATE( L1error(1:nVar) )
ALLOCATE( Wexact(1:nVar,1:nElemsX,1:nElemsY+1)) !*NB:+1 in Y direction
 
Wexact(:,:,:) = 0.
L1error(:) = 0. 
WtempCons = 0.
WtempPrim = 0.
DO jj = 1,nElemsY+1 !*NB:+1 in Y direction
  DO ii = 1,nElemsX
    DO iGP=1,nGPs
      DO jGP=1,nGPs
        CALL ExactFunction(InitialCondition, tEnd, MeshGP_Y(:,ii,jj,iGP,jGP), WtempCons(1:nVar,iGP,jGP))
        CALL ConsToPrim(WtempCons(1:nVar,iGP,jGP),WtempPrim(1:nVar,iGP,jGP))
        Wexact(1:nVar,ii,jj) = Wexact(1:nVar,ii,jj) + WeightsGP(iGP,jGP)* WtempPrim(1:nVar,iGP,jGP)
      END DO
    END DO
    L1error(:) = L1error(:) + ABS( W_Y(:,ii,jj) - Wexact(:,ii,jj) ) 
  ENDDO
ENDDO

L1error(:) = L1error(:) * MESH_DX(1) * MESH_DX(2)

PRINT*, "Error W_Y", L1error

FNAME = 'W_Y_ErrorL1_XXXX_XXXX.dat'
WRITE(FNAME(9+4:12+4),FMT="(I4.4)") nElemsX
WRITE(FNAME(14+4:17+4),FMT="(I4.4)") nElemsY

OPEN(666,FILE=TRIM(FNAME))
WRITE(666,*) nElemsX, nElemsY, L1error(1:nVar)
#ifdef ACTIVEFLUX
#ifdef IMEX
WRITE(666,*) HyperbolicityLoss
#endif
#endif

DEALLOCATE( L1error )
DEALLOCATE( Wexact )

#endif

#ifdef PATANKAR
averageJacobi = REAL(SUM(JacobiIterations(1:JacobiCounter)))/JacobiCounter
std_deviation = SQRT(SUM((REAL(JacobiIterations(1:JacobiCounter))-averageJacobi)**2)/(JacobiCounter-1))
WRITE(666,*) 'Jacobi average iterations'
WRITE(666,*) averageJacobi
WRITE(666,*) 'Jacobi standard deviation'
WRITE(666,*) std_deviation
#endif

CLOSE(666)

#ifndef PRIMITIVEONLY
open(unit=7,file='den251')
i=0
do k=1,nElemsY
  do j=1,nElemsX
      i=i+1
      write(7,*) i, U(1,j,k)
  enddo
enddo
#endif

#ifdef CENTEREDPRIMITIVE
ALLOCATE( L1error(1:nVar) )
ALLOCATE( Wexact(1:nVar,1:nElemsX,1:nElemsY)) 

Wexact(:,:,:) = 0.
L1error(:) = 0. 
WtempCons = 0.
WtempPrim = 0.
DO jj = 1,nElemsY
  DO ii = 1,nElemsX
    DO iGP=1,nGPs
      DO jGP=1,nGPs
        CALL ExactFunction(InitialCondition, tEnd, MeshGP(:,ii,jj,iGP,jGP), WtempCons(1:nVar,iGP,jGP))
        CALL ConsToPrim(WtempCons(1:nVar,iGP,jGP),WtempPrim(1:nVar,iGP,jGP))
        Wexact(1:nVar,ii,jj) = Wexact(1:nVar,ii,jj) + WeightsGP(iGP,jGP)* WtempPrim(1:nVar,iGP,jGP)
      END DO
    END DO
    L1error(:) = L1error(:) + ABS( WC(:,ii,jj) - Wexact(:,ii,jj) ) 
  ENDDO
ENDDO

L1error(:) = L1error(:) * MESH_DX(1) * MESH_DX(2)

PRINT*, "Error WC", L1error

FNAME = 'WC_ErrorL1_XXXX_XXXX.dat'
WRITE(FNAME(9+3:12+3),FMT="(I4.4)") nElemsX
WRITE(FNAME(14+3:17+3),FMT="(I4.4)") nElemsY

OPEN(666,FILE=TRIM(FNAME))
WRITE(666,*) nElemsX, nElemsY, L1error(1:nVar)
#endif
#ifdef CENTEREDPRIMITIVE
#ifdef IMEX
WRITE(666,*) HyperbolicityLoss
#endif
#endif

PRINT*, "Computational time", computationalTime
PRINT*, "--------------------------"

END SUBROUTINE

END MODULE MOD_PostProcessing
