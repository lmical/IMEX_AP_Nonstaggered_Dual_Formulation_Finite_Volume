!===============================================================================!
MODULE MOD_Reconstruction
!===============================================================================!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE ReconstructionX
  MODULE PROCEDURE ReconstructionX
END INTERFACE

INTERFACE ReconstructionY
  MODULE PROCEDURE ReconstructionY
END INTERFACE

INTERFACE ReconstructionFixX
  MODULE PROCEDURE ReconstructionFixX
END INTERFACE

INTERFACE ReconstructionFixY
  MODULE PROCEDURE ReconstructionFixY
END INTERFACE

INTERFACE MUSCL 
  MODULE PROCEDURE MUSCL 
END INTERFACE

INTERFACE kMUSCL 
  MODULE PROCEDURE kMUSCL 
END INTERFACE

INTERFACE COMUSCL 
  MODULE PROCEDURE COMUSCL 
END INTERFACE

INTERFACE VLMUSCL 
  MODULE PROCEDURE VLMUSCL 
END INTERFACE

INTERFACE M2MUSCL 
  MODULE PROCEDURE M2MUSCL 
END INTERFACE

INTERFACE VAMUSCL 
  MODULE PROCEDURE VAMUSCL 
END INTERFACE

INTERFACE SBMMUSCL 
  MODULE PROCEDURE SBMMUSCL 
END INTERFACE

INTERFACE THETAMUSCL 
  MODULE PROCEDURE THETAMUSCL 
END INTERFACE

INTERFACE WENO3_SecondSweep 
  MODULE PROCEDURE WENO3_SecondSweep
END INTERFACE

INTERFACE WENO5_SecondSweep 
  MODULE PROCEDURE WENO5_SecondSweep
END INTERFACE

INTERFACE WENO7_SecondSweep 
  MODULE PROCEDURE WENO7_SecondSweep
END INTERFACE

INTERFACE WENO_In_Quadrature_Points 
  MODULE PROCEDURE WENO_In_Quadrature_Points
END INTERFACE

INTERFACE First_Derivative_Central_Order2 
  MODULE PROCEDURE First_Derivative_Central_Order2
END INTERFACE

INTERFACE Second_Derivative_Central_Order2 
  MODULE PROCEDURE Second_Derivative_Central_Order2
END INTERFACE

INTERFACE Second_Mixed_Derivative_Central_Order2 
  MODULE PROCEDURE Second_Mixed_Derivative_Central_Order2
END INTERFACE

INTERFACE Laplacian_Central_Order2 
  MODULE PROCEDURE Laplacian_Central_Order2
END INTERFACE

!*Unified
INTERFACE ReconstructionX_INPUT_PRIMITIVE
  MODULE PROCEDURE ReconstructionX_INPUT_PRIMITIVE
END INTERFACE

INTERFACE ReconstructionY_INPUT_PRIMITIVE
  MODULE PROCEDURE ReconstructionY_INPUT_PRIMITIVE
END INTERFACE

INTERFACE ReconstructionX_INPUT_CONSERVED
  MODULE PROCEDURE ReconstructionX_INPUT_CONSERVED
END INTERFACE

INTERFACE ReconstructionY_INPUT_CONSERVED
  MODULE PROCEDURE ReconstructionY_INPUT_CONSERVED
END INTERFACE

INTERFACE ReconstructionFixX_INPUT_PRIMITIVE
  MODULE PROCEDURE ReconstructionFixX_INPUT_PRIMITIVE
END INTERFACE

INTERFACE ReconstructionFixY_INPUT_PRIMITIVE
  MODULE PROCEDURE ReconstructionFixY_INPUT_PRIMITIVE
END INTERFACE

INTERFACE ReconstructionFixX_INPUT_CONSERVED
  MODULE PROCEDURE ReconstructionFixX_INPUT_CONSERVED
END INTERFACE

INTERFACE ReconstructionFixY_INPUT_CONSERVED
  MODULE PROCEDURE ReconstructionFixY_INPUT_CONSERVED
END INTERFACE


INTERFACE MINMOD 
  MODULE PROCEDURE MINMOD
END INTERFACE

INTERFACE MINMOD_3_INPUTS 
  MODULE PROCEDURE MINMOD_3_INPUTS
END INTERFACE

#ifdef ACTIVEFLUX
INTERFACE FromCellAveragesToPointValues 
  MODULE PROCEDURE FromCellAveragesToPointValues
END INTERFACE
#endif
!-------------------------------------------------------------------------------!
PUBLIC :: ReconstructionX
PUBLIC :: ReconstructionY
PUBLIC :: ReconstructionFixX
PUBLIC :: ReconstructionFixY
PUBLIC :: MUSCL
PUBLIC :: WENO3_SecondSweep 
PUBLIC :: WENO5_SecondSweep 
PUBLIC :: WENO7_SecondSweep 
PUBLIC :: WENO_In_Quadrature_Points
PUBLIC :: First_Derivative_Central_Order2
PUBLIC :: Second_Derivative_Central_Order2
PUBLIC :: Second_Mixed_Derivative_Central_Order2
PUBLIC :: Laplacian_Central_Order2
!*Unified
PUBLIC :: ReconstructionX_INPUT_PRIMITIVE
PUBLIC :: ReconstructionY_INPUT_PRIMITIVE
PUBLIC :: ReconstructionX_INPUT_CONSERVED
PUBLIC :: ReconstructionY_INPUT_CONSERVED
PUBLIC :: ReconstructionFixX_INPUT_PRIMITIVE
PUBLIC :: ReconstructionFixY_INPUT_PRIMITIVE
PUBLIC :: ReconstructionFixX_INPUT_CONSERVED
PUBLIC :: ReconstructionFixY_INPUT_CONSERVED
PUBLIC :: MINMOD 
PUBLIC :: MINMOD_3_INPUTS
PUBLIC :: kMUSCL
PUBLIC :: COMUSCL
PUBLIC :: VLMUSCL
PUBLIC :: M2MUSCL
PUBLIC :: VAMUSCL
PUBLIC :: SBMMUSCL
PUBLIC :: THETAMUSCL

#ifdef ACTIVEFLUX
PUBLIC :: FromCellAveragesToPointValues 

!*AF Postprocessing
PUBLIC :: AF_MUSCL
PUBLIC :: AF_MUSCL_Primitive_Based
PUBLIC :: SIDED_AF_MUSCL
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
SUBROUTINE ReconstructionX()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: LLX
USE MOD_FiniteVolume2D_vars,ONLY: RRX
USE MOD_FiniteVolume2D_vars,ONLY: LLY
USE MOD_FiniteVolume2D_vars,ONLY: RRY
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
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructedVariable
USE MOD_FiniteVolume2D_vars,ONLY: Tau_Standard_Limiting
USE MOD_FiniteVolume2D_vars,ONLY: Theta_Standard_Limiting
USE MOD_FiniteVolume2D_vars,ONLY: Tau_Overcompressive
USE MOD_FiniteVolume2D_vars,ONLY: Theta_Overcompressive
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj, kk, iGP
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: W(1:nVar,-nGhosts:nGhosts)
!-------------------------------------------------------------------------------!
DO jj=1,nElemsY
  DO ii=0,nElemsX+1

    SELECT CASE (ReconstructedVariable)
      CASE(0)
        W(1:nVar,-nGhosts:nGhosts)=U(1:nVar,-nGhosts+ii:ii+nGhosts,jj)
      CASE(1)
        DO kk=-nGhosts,nGhosts
          W(1:nVar,kk)=MATMUL(LLX(:,:,ii,jj),U(1:nVar,ii+kk,jj))
        END DO
      CASE DEFAULT
        ErrorMessage = "Error in ReconstructionX. ReconstructedVariable not available"
        WRITE(*,*) ErrorMessage
        STOP
    END SELECT


    SELECT CASE (ABS(Reconstruction))
      CASE(1,10)
            WM(1:nVar,nGPs,ii,jj) = W(1:nVar,0)
            WP(1:nVar,nGPs,ii,jj) = W(1:nVar,0)
      CASE(2,20) !*MUSCL
            CALL MUSCL(&
                      W(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
      CASE(21) !*k-MUSCL
            CALL kMUSCL(&
                      W(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
      CASE(22) !*MUSCL_CO
            CALL COMUSCL(&
                      W(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
      CASE(23) !*MUSCL_VL
            CALL VLMUSCL(&
                      W(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
      CASE(24) !*MUSCL_M
            CALL M2MUSCL(&
                      W(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
      CASE(25) !*MUSCL_VA
            CALL VAMUSCL(&
                      W(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
      CASE(26) !*MUSCL_SBM
            CALL SBMMUSCL(&
                      W(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
      CASE(27) !*MUSCL_THETA
            CALL THETAMUSCL(&
                      W(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
      CASE(28) !*MUSCL_SBM standard
            CALL SBMMUSCL(&
                      W(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1),tau_input=Tau_Standard_Limiting,theta_input=Theta_Standard_Limiting)
      CASE(29) !*MUSCL_SBM overcompressive
            CALL SBMMUSCL(&
                      W(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1),tau_input=Tau_Overcompressive,theta_input=Theta_Overcompressive)
      CASE(3,4,5,7)
            IF (.NOT. Ind(1,ii,jj)) THEN
              !*NB: Here I'm entering with conserved variables
              CALL WENO_XDIR(&
                        U(1:nVar,-nGhosts+ii:ii+nGhosts,-nGhosts+jj:jj+nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(1),&
                        Reconstruction,&
                        RRX(1:nVar,1:nVar,ii,jj),&
                        LLX(1:nVar,1:nVar,ii,jj),&
                        RRY(1:nVar,1:nVar,ii,jj),&
                        LLY(1:nVar,1:nVar,ii,jj))
            END IF
      CASE DEFAULT
        ErrorMessage = "Reconstruction not implemented"
        WRITE(*,*) ErrorMessage
        STOP
    END SELECT

    !*We transform back
    IF ((ABS(Reconstruction) .LT. 3) .OR. (ABS(Reconstruction) .EQ. 10) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 20) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 21) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 22) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 23) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 24) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 25) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 26) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 27) &                                  
                                  &  .OR. (ABS(Reconstruction) .EQ. 28) &                                  
                                  &  .OR. (ABS(Reconstruction) .EQ. 29) ) THEN 
      SELECT CASE (ReconstructedVariable)
        CASE(0)

        CASE(1)
          DO iGP=1,nGPs
            WM(1:nVar,iGP,ii,jj)=MATMUL(RRX(:,:,ii,jj),WM(1:nVar,iGP,ii,jj))
            WP(1:nVar,iGP,ii,jj)=MATMUL(RRX(:,:,ii,jj),WP(1:nVar,iGP,ii,jj))
          END DO
        CASE DEFAULT
          ErrorMessage = "Error in ReconstructionX. ReconstructedVariable not available"
          WRITE(*,*) ErrorMessage
          STOP
      END SELECT
    END IF

  END DO
END DO
!-------------------------------------------------------------------------------!
END SUBROUTINE ReconstructionX
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ReconstructionY()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: LLX
USE MOD_FiniteVolume2D_vars,ONLY: RRX
USE MOD_FiniteVolume2D_vars,ONLY: LLY
USE MOD_FiniteVolume2D_vars,ONLY: RRY
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructedVariable
USE MOD_FiniteVolume2D_vars,ONLY: Tau_Standard_Limiting
USE MOD_FiniteVolume2D_vars,ONLY: Theta_Standard_Limiting
USE MOD_FiniteVolume2D_vars,ONLY: Tau_Overcompressive
USE MOD_FiniteVolume2D_vars,ONLY: Theta_Overcompressive
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj, kk, iGP
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: W(1:nVar,-nGhosts:nGhosts)
!-------------------------------------------------------------------------------!

DO jj=0,nElemsY+1
  DO ii=1,nElemsX

    SELECT CASE (ReconstructedVariable)
      CASE(0)
        W(1:nVar,-nGhosts:nGhosts)=U(1:nVar,ii,-nGhosts+jj:jj+nGhosts)
      CASE(1)
        DO kk=-nGhosts,nGhosts
          W(1:nVar,kk)=MATMUL(LLY(:,:,ii,jj),U(1:nVar,ii,jj+kk))
        END DO
      CASE DEFAULT
        ErrorMessage = "Error in ReconstructionY. ReconstructedVariable not available"
        WRITE(*,*) ErrorMessage
        STOP
    END SELECT



    SELECT CASE (ABS(Reconstruction))
      CASE(1,10)
            WM(1:nVar,nGPs,ii,jj) = W(1:nVar,0)
            WP(1:nVar,nGPs,ii,jj) = W(1:nVar,0)
      CASE(2,20) !*MUSCL
      CALL MUSCL(&
                W(1:nVar,-nGhosts:nGhosts),&
                WM(1:nVar,1:nGPs,ii,jj),&
                WP(1:nVar,1:nGPs,ii,jj),&
                MESH_DX(2))
      CASE(21) !*k-MUSCL
        CALL kMUSCL(&
                  W(1:nVar,-nGhosts:nGhosts),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(2))
      CASE(22) !*MUSCL_CO
        CALL COMUSCL(&
                  W(1:nVar,-nGhosts:nGhosts),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(2))
      CASE(23) !*MUSCL_VL
        CALL VLMUSCL(&
                  W(1:nVar,-nGhosts:nGhosts),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(2))
      CASE(24) !*MUSCL_M
        CALL M2MUSCL(&
                  W(1:nVar,-nGhosts:nGhosts),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(2))
      CASE(25) !*MUSCL_VA
        CALL VAMUSCL(&
                  W(1:nVar,-nGhosts:nGhosts),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(2))
      CASE(26) !*MUSCL_SBM
        CALL SBMMUSCL(&
                  W(1:nVar,-nGhosts:nGhosts),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(2))
      CASE(27) !*MUSCL_THETA
        CALL THETAMUSCL(&
                  W(1:nVar,-nGhosts:nGhosts),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(2))
      CASE(28) !*MUSCL_SBM standard
            CALL SBMMUSCL(&
                  W(1:nVar,-nGhosts:nGhosts),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(2),tau_input=Tau_Standard_Limiting,theta_input=Theta_Standard_Limiting)
      CASE(29) !*MUSCL_SBM overcompressive
            CALL SBMMUSCL(&
                  W(1:nVar,-nGhosts:nGhosts),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(2),tau_input=Tau_Overcompressive,theta_input=Theta_Overcompressive)
      CASE(3,4,5,7)
        IF (.NOT. Ind(2,ii,jj)) THEN
          !*NB: Here I'm entering with conserved variables
          CALL WENO_YDIR(&
                    U(1:nVar,-nGhosts+ii:ii+nGhosts,-nGhosts+jj:jj+nGhosts),&
                    WM(1:nVar,1:nGPs,ii,jj),&
                    WP(1:nVar,1:nGPs,ii,jj),&
                    MESH_DX(2),&
                    Reconstruction,&
                    RRX(1:nVar,1:nVar,ii,jj),&
                    LLX(1:nVar,1:nVar,ii,jj),&
                    RRY(1:nVar,1:nVar,ii,jj),&
                    LLY(1:nVar,1:nVar,ii,jj))
        END IF
      CASE DEFAULT
        ErrorMessage = "Reconstruction not implemented"
        WRITE(*,*) ErrorMessage
        STOP
    END SELECT


    !*We transform back
    IF ((ABS(Reconstruction) .LT. 3) .OR. (ABS(Reconstruction) .EQ. 10) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 20) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 21) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 22) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 23) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 24) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 25) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 26) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 27) &                                  
                                  &  .OR. (ABS(Reconstruction) .EQ. 28) &                                  
                                  &  .OR. (ABS(Reconstruction) .EQ. 29) ) THEN 
      SELECT CASE (ReconstructedVariable)
        CASE(0)

        CASE(1)
          DO iGP=1,nGPs
            WM(1:nVar,iGP,ii,jj)=MATMUL(RRY(:,:,ii,jj),WM(1:nVar,iGP,ii,jj))
            WP(1:nVar,iGP,ii,jj)=MATMUL(RRY(:,:,ii,jj),WP(1:nVar,iGP,ii,jj))
          END DO
        CASE DEFAULT
          ErrorMessage = "Error in ReconstructionY. ReconstructedVariable not available"
          WRITE(*,*) ErrorMessage
          STOP
      END SELECT
    END IF


  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE ReconstructionY
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ReconstructionFixX()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: LLX
USE MOD_FiniteVolume2D_vars,ONLY: RRX
USE MOD_FiniteVolume2D_vars,ONLY: LLY
USE MOD_FiniteVolume2D_vars,ONLY: RRY
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructionFix
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructedVariable
USE MOD_FiniteVolume2D_vars,ONLY: Tau_Standard_Limiting
USE MOD_FiniteVolume2D_vars,ONLY: Theta_Standard_Limiting
USE MOD_FiniteVolume2D_vars,ONLY: Tau_Overcompressive
USE MOD_FiniteVolume2D_vars,ONLY: Theta_Overcompressive
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj, kk, iGP
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: W(1:nVar,-nGhosts:nGhosts)
!-------------------------------------------------------------------------------!

SELECT CASE (ABS(Reconstruction))
  CASE(1,10,2,20,21,22,23,24,25,26,27,28,29)
    RETURN
  CASE DEFAULT
    !*Nothing
END SELECT  


DO jj=1,nElemsY
  DO ii=0,nElemsX+1

    IF (Ind(1,ii,jj)) THEN

      SELECT CASE (ReconstructedVariable)
        CASE(0)
          W(1:nVar,-nGhosts:nGhosts)=U(1:nVar,-nGhosts+ii:ii+nGhosts,jj)
        CASE(1)
          DO kk=-nGhosts,nGhosts
            W(1:nVar,kk)=MATMUL(LLX(:,:,ii,jj),U(1:nVar,ii+kk,jj))
          END DO
        CASE DEFAULT
          ErrorMessage = "Error in ReconstructionFixX. ReconstructedVariable not available"
          WRITE(*,*) ErrorMessage
          STOP
      END SELECT



      SELECT CASE (ABS(ReconstructionFix))
        CASE(1,10)
          DO iGP=1,nGPs
              WM(1:nVar,iGP,ii,jj) = W(1:nVar,0)
              WP(1:nVar,iGP,ii,jj) = W(1:nVar,0)
          END DO
        CASE(2,20)
          CALL MUSCL(&
                    W(1:nVar,-nGhosts:nGhosts),&
                    WM(1:nVar,1:nGPs,ii,jj),&
                    WP(1:nVar,1:nGPs,ii,jj),&
                    MESH_DX(1))
        CASE(21) !*k-MUSCL
          CALL kMUSCL(&
                    W(1:nVar,-nGhosts:nGhosts),&
                    WM(1:nVar,1:nGPs,ii,jj),&
                    WP(1:nVar,1:nGPs,ii,jj),&
                    MESH_DX(1))
        CASE(22) !*MUSCL_CO
          CALL COMUSCL(&
                    W(1:nVar,-nGhosts:nGhosts),&
                    WM(1:nVar,1:nGPs,ii,jj),&
                    WP(1:nVar,1:nGPs,ii,jj),&
                    MESH_DX(1))
        CASE(23) !*MUSCL_VL
          CALL VLMUSCL(&
                    W(1:nVar,-nGhosts:nGhosts),&
                    WM(1:nVar,1:nGPs,ii,jj),&
                    WP(1:nVar,1:nGPs,ii,jj),&
                    MESH_DX(1))
        CASE(24) !*MUSCL_M
          CALL M2MUSCL(&
                    W(1:nVar,-nGhosts:nGhosts),&
                    WM(1:nVar,1:nGPs,ii,jj),&
                    WP(1:nVar,1:nGPs,ii,jj),&
                    MESH_DX(1))
        CASE(25) !*MUSCL_VA
          CALL VAMUSCL(&
                    W(1:nVar,-nGhosts:nGhosts),&
                    WM(1:nVar,1:nGPs,ii,jj),&
                    WP(1:nVar,1:nGPs,ii,jj),&
                    MESH_DX(1))
        CASE(26) !*MUSCL_SBM
          CALL SBMMUSCL(&
                    W(1:nVar,-nGhosts:nGhosts),&
                    WM(1:nVar,1:nGPs,ii,jj),&
                    WP(1:nVar,1:nGPs,ii,jj),&
                    MESH_DX(1))
        CASE(27) !*MUSCL_THETA
          CALL THETAMUSCL(&
                    W(1:nVar,-nGhosts:nGhosts),&
                    WM(1:nVar,1:nGPs,ii,jj),&
                    WP(1:nVar,1:nGPs,ii,jj),&
                    MESH_DX(1))
      CASE(28) !*MUSCL_SBM standard
            CALL SBMMUSCL(&
                      W(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1),tau_input=Tau_Standard_Limiting,theta_input=Theta_Standard_Limiting)
      CASE(29) !*MUSCL_SBM overcompressive
            CALL SBMMUSCL(&
                      W(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1),tau_input=Tau_Overcompressive,theta_input=Theta_Overcompressive)
        CASE(3,4,5,7)
          !*NB: Here, I'm entering with conseved variables
          CALL WENO_XDIR(&
                    U(1:nVar,-nGhosts+ii:ii+nGhosts,-nGhosts+jj:jj+nGhosts),&
                    WM(1:nVar,1:nGPs,ii,jj),&
                    WP(1:nVar,1:nGPs,ii,jj),&
                    MESH_DX(1),&
                    ReconstructionFix,&
                    RRX(1:nVar,1:nVar,ii,jj),&
                    LLX(1:nVar,1:nVar,ii,jj),&
                    RRY(1:nVar,1:nVar,ii,jj),&
                    LLY(1:nVar,1:nVar,ii,jj))
        CASE DEFAULT
          ErrorMessage = "ReconstructionFixX not implemented"
          WRITE(*,*) ErrorMessage
          STOP
      END SELECT



    !*We transform back
    IF ((ABS(Reconstruction) .LT. 3) .OR. (ABS(Reconstruction) .EQ. 10) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 20) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 21) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 22) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 23) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 24) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 25) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 26) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 27) &                                  
                                  &  .OR. (ABS(Reconstruction) .EQ. 28) &                                  
                                  &  .OR. (ABS(Reconstruction) .EQ. 29) ) THEN 
        SELECT CASE (ReconstructedVariable)
          CASE(0)

          CASE(1)
            DO iGP=1,nGPs
              WM(1:nVar,iGP,ii,jj)=MATMUL(RRX(:,:,ii,jj),WM(1:nVar,iGP,ii,jj))
              WP(1:nVar,iGP,ii,jj)=MATMUL(RRX(:,:,ii,jj),WP(1:nVar,iGP,ii,jj))
            END DO
          CASE DEFAULT
            ErrorMessage = "Error in ReconstructionFixX. ReconstructedVariable not available"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT
      END IF




    END IF

  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE ReconstructionFixX
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ReconstructionFixY()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: LLX
USE MOD_FiniteVolume2D_vars,ONLY: RRX
USE MOD_FiniteVolume2D_vars,ONLY: LLY
USE MOD_FiniteVolume2D_vars,ONLY: RRY
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructionFix
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructedVariable
USE MOD_FiniteVolume2D_vars,ONLY: Tau_Standard_Limiting
USE MOD_FiniteVolume2D_vars,ONLY: Theta_Standard_Limiting
USE MOD_FiniteVolume2D_vars,ONLY: Tau_Overcompressive
USE MOD_FiniteVolume2D_vars,ONLY: Theta_Overcompressive
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj, kk, iGP
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: W(1:nVar,-nGhosts:nGhosts)
!-------------------------------------------------------------------------------!

SELECT CASE (ABS(Reconstruction))
  CASE(1,10,2,20,21,22,23,24,25,26,27,28,29)
    RETURN
  CASE DEFAULT
    !*Nothing
END SELECT  


DO jj=0,nElemsY+1
  DO ii=1,nElemsX
    IF (Ind(2,ii,jj)) THEN

      SELECT CASE (ReconstructedVariable)
        CASE(0)
          W(1:nVar,-nGhosts:nGhosts)=U(1:nVar,ii,-nGhosts+jj:jj+nGhosts)
        CASE(1)
          DO kk=-nGhosts,nGhosts
            W(1:nVar,kk)=MATMUL(LLY(:,:,ii,jj),U(1:nVar,ii,jj+kk))
          END DO
        CASE DEFAULT
          ErrorMessage = "Error in ReconstructionFixY. ReconstructedVariable not available"
          WRITE(*,*) ErrorMessage
          STOP
      END SELECT


      SELECT CASE (ABS(ReconstructionFix))
        CASE(1,10)
          DO iGP=1,nGPs
              WM(1:nVar,iGP,ii,jj) = W(1:nVar,0)
              WP(1:nVar,iGP,ii,jj) = W(1:nVar,0)
          END DO
        CASE(2,20)
          CALL MUSCL(&
                    W(1:nVar,-nGhosts:+nGhosts),&
                    WM(1:nVar,1:nGPs,ii,jj),&
                    WP(1:nVar,1:nGPs,ii,jj),&
                    MESH_DX(2))
        CASE(21) !*k-MUSCL
          CALL kMUSCL(&
                    W(1:nVar,-nGhosts:+nGhosts),&
                    WM(1:nVar,1:nGPs,ii,jj),&
                    WP(1:nVar,1:nGPs,ii,jj),&
                    MESH_DX(2))
        CASE(22) !*MUSCL_CO
          CALL COMUSCL(&
                    W(1:nVar,-nGhosts:+nGhosts),&
                    WM(1:nVar,1:nGPs,ii,jj),&
                    WP(1:nVar,1:nGPs,ii,jj),&
                    MESH_DX(2))
        CASE(23) !*MUSCL_VL
          CALL VLMUSCL(&
                    W(1:nVar,-nGhosts:+nGhosts),&
                    WM(1:nVar,1:nGPs,ii,jj),&
                    WP(1:nVar,1:nGPs,ii,jj),&
                    MESH_DX(2))
        CASE(24) !*MUSCL_M
          CALL M2MUSCL(&
                    W(1:nVar,-nGhosts:+nGhosts),&
                    WM(1:nVar,1:nGPs,ii,jj),&
                    WP(1:nVar,1:nGPs,ii,jj),&
                    MESH_DX(2))
        CASE(25) !*MUSCL_VA
          CALL VAMUSCL(&
                    W(1:nVar,-nGhosts:+nGhosts),&
                    WM(1:nVar,1:nGPs,ii,jj),&
                    WP(1:nVar,1:nGPs,ii,jj),&
                    MESH_DX(2))
        CASE(26) !*MUSCL_SBM
          CALL SBMMUSCL(&
                    W(1:nVar,-nGhosts:+nGhosts),&
                    WM(1:nVar,1:nGPs,ii,jj),&
                    WP(1:nVar,1:nGPs,ii,jj),&
                    MESH_DX(2))
        CASE(27) !*MUSCL_THETA
          CALL THETAMUSCL(&
                    W(1:nVar,-nGhosts:+nGhosts),&
                    WM(1:nVar,1:nGPs,ii,jj),&
                    WP(1:nVar,1:nGPs,ii,jj),&
                    MESH_DX(2))
      CASE(28) !*MUSCL_SBM standard
            CALL SBMMUSCL(&
                  W(1:nVar,-nGhosts:nGhosts),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(2),tau_input=Tau_Standard_Limiting,theta_input=Theta_Standard_Limiting)
      CASE(29) !*MUSCL_SBM overcompressive
            CALL SBMMUSCL(&
                  W(1:nVar,-nGhosts:nGhosts),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(2),tau_input=Tau_Overcompressive,theta_input=Theta_Overcompressive)
        CASE(3,4,5,7)
          !*NB: Here, I'm entering with conseved variables
          CALL WENO_YDIR(&
                    U(1:nVar,-nGhosts+ii:ii+nGhosts,-nGhosts+jj:jj+nGhosts),&
                    WM(1:nVar,1:nGPs,ii,jj),&
                    WP(1:nVar,1:nGPs,ii,jj),&
                    MESH_DX(2),&
                    ReconstructionFix,&
                    RRX(1:nVar,1:nVar,ii,jj),&
                    LLX(1:nVar,1:nVar,ii,jj),&
                    RRY(1:nVar,1:nVar,ii,jj),&
                    LLY(1:nVar,1:nVar,ii,jj))
        CASE DEFAULT
          ErrorMessage = "ReconstructionFix not implemented"
          WRITE(*,*) ErrorMessage
          STOP
      END SELECT


    !*We transform back
    IF ((ABS(Reconstruction) .LT. 3) .OR. (ABS(Reconstruction) .EQ. 10) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 20) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 21) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 22) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 23) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 24) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 25) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 26) &
                                  &  .OR. (ABS(Reconstruction) .EQ. 27) &                                  
                                  &  .OR. (ABS(Reconstruction) .EQ. 28) &                                  
                                  &  .OR. (ABS(Reconstruction) .EQ. 29) ) THEN 
        SELECT CASE (ReconstructedVariable)
          CASE(0)

          CASE(1)
            DO iGP=1,nGPs
              WM(1:nVar,iGP,ii,jj)=MATMUL(RRY(:,:,ii,jj),WM(1:nVar,iGP,ii,jj))
              WP(1:nVar,iGP,ii,jj)=MATMUL(RRY(:,:,ii,jj),WP(1:nVar,iGP,ii,jj))
            END DO
          CASE DEFAULT
            ErrorMessage = "Error in ReconstructionY. ReconstructedVariable not available"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT
      END IF




    END IF
  END DO
END DO


!-------------------------------------------------------------------------------!
END SUBROUTINE ReconstructionFixY
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE MUSCL(Q,WM,WP,dx)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: dx
REAL,INTENT(IN)  :: Q(1:nVar,-nGhosts:nGhosts)
REAL,INTENT(OUT) :: WM(1:nVar,1:nGPs)
REAL,INTENT(OUT) :: WP(1:nVar,1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: sm, sp, slope
INTEGER          :: iVar, iGP
!-------------------------------------------------------------------------------!

DO iGP=1,nGPs
  DO iVar = 1, nVar
    sp    = (Q(iVar,+1) - Q(iVar,+0))/dx
    sm    = (Q(iVar,+0) - Q(iVar,-1))/dx
    slope = MINMOD(sm,sp)

    WM(iVar,iGP) = Q(iVar,0) - 0.5*slope*dx
    WP(iVar,iGP) = Q(iVar,0) + 0.5*slope*dx
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE MUSCL
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE kMUSCL(Q,WM,WP,dx)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts

!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: dx
REAL,INTENT(IN)  :: Q(1:nVar,-nGhosts:nGhosts)
REAL,INTENT(OUT) :: WM(1:nVar,1:nGPs)
REAL,INTENT(OUT) :: WP(1:nVar,1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: slope, a, b, s1, s2
REAL, PARAMETER  :: k=1.5 !* 1<=k<=2
INTEGER          :: iVar, iGP
!-------------------------------------------------------------------------------!

DO iGP=1,nGPs
  DO iVar = 1, nVar
    a    = Q(iVar,+1) - Q(iVar,+0)
    b    = Q(iVar,+0) - Q(iVar,-1)

    s1=ABS( MINMOD( k*a , b   ) )
    s2=ABS( MINMOD( a   , k*b ) )

    !*SIGN(A,B) returns the value of A with the sign of B.
    !*NB: SIGN should work for reals
    slope = SIGN(1.,a)/dx*MAX( s1, s2 )
      
    WM(iVar,iGP) = Q(iVar,0) - 0.5*slope*dx
    WP(iVar,iGP) = Q(iVar,0) + 0.5*slope*dx
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE kMUSCL
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE COMUSCL(Q,WM,WP,dx)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: dx
REAL,INTENT(IN)  :: Q(1:nVar,-nGhosts:nGhosts)
REAL,INTENT(OUT) :: WM(1:nVar,1:nGPs)
REAL,INTENT(OUT) :: WP(1:nVar,1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: slope, a, b, s1
REAL, PARAMETER  :: k=1. !* 1<=k<=2 !*NB: In my opinion it only works with k=1
INTEGER          :: iVar, iGP
!-------------------------------------------------------------------------------!

DO iGP=1,nGPs
  DO iVar = 1, nVar
    a    = Q(iVar,+1) - Q(iVar,+0)
    b    = Q(iVar,+0) - Q(iVar,-1)

    s1=MINMOD( k*a , b )
    slope = s1/dx
      
    WM(iVar,iGP) = Q(iVar,0) - 0.5*slope*dx
    WP(iVar,iGP) = Q(iVar,0) + 0.5*slope*dx
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE COMUSCL
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE VLMUSCL(Q,WM,WP,dx)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: dx
REAL,INTENT(IN)  :: Q(1:nVar,-nGhosts:nGhosts)
REAL,INTENT(OUT) :: WM(1:nVar,1:nGPs)
REAL,INTENT(OUT) :: WP(1:nVar,1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: slope, a, b
INTEGER          :: iVar, iGP
!-------------------------------------------------------------------------------!

DO iGP=1,nGPs
  DO iVar = 1, nVar
    a    = Q(iVar,+1) - Q(iVar,+0)
    b    = Q(iVar,+0) - Q(iVar,-1)

    slope = 0.

    IF ( (a*b) .GT. 0. ) THEN
      slope = 2.*a*b/(a+b)/dx
    ELSE
      slope = 0.
    END IF

    WM(iVar,iGP) = Q(iVar,0) - 0.5*slope*dx
    WP(iVar,iGP) = Q(iVar,0) + 0.5*slope*dx
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE VLMUSCL
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE M2MUSCL(Q,WM,WP,dx)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: dx
REAL,INTENT(IN)  :: Q(1:nVar,-nGhosts:nGhosts)
REAL,INTENT(OUT) :: WM(1:nVar,1:nGPs)
REAL,INTENT(OUT) :: WP(1:nVar,1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: slope, a, b, s1, s2
INTEGER          :: iVar, iGP
!-------------------------------------------------------------------------------!

DO iGP=1,nGPs
  DO iVar = 1, nVar
    a    = Q(iVar,+1) - Q(iVar,+0)
    b    = Q(iVar,+0) - Q(iVar,-1)

    s1=(a+b)/2.
    s2=2.*MINMOD(a,b)

    slope = MINMOD(s1,s2)/dx

    WM(iVar,iGP) = Q(iVar,0) - 0.5*slope*dx
    WP(iVar,iGP) = Q(iVar,0) + 0.5*slope*dx
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE M2MUSCL
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE VAMUSCL(Q,WM,WP,dx)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: dx
REAL,INTENT(IN)  :: Q(1:nVar,-nGhosts:nGhosts)
REAL,INTENT(OUT) :: WM(1:nVar,1:nGPs)
REAL,INTENT(OUT) :: WP(1:nVar,1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: slope, a, b, c
INTEGER          :: iVar, iGP
!-------------------------------------------------------------------------------!

DO iGP=1,nGPs
  DO iVar = 1, nVar
    a    = Q(iVar,+1) - Q(iVar,+0)
    b    = Q(iVar,+0) - Q(iVar,-1)

    PRINT*, "MUSCL_VA coded but I do not really know what the parameter c is"
    PRINT*, "Check before using it!"
    STOP

    slope = ( a*b + c**2 ) * (a+b) / (a**2 + b**2 + 2.*c**2 ) / dx

    WM(iVar,iGP) = Q(iVar,0) - 0.5*slope*dx
    WP(iVar,iGP) = Q(iVar,0) + 0.5*slope*dx
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE VAMUSCL
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE SBMMUSCL(Q,WM,WP,dx,slope_output,tau_input,theta_input)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)              :: dx
REAL,INTENT(IN)              :: Q(1:nVar,-nGhosts:nGhosts)
REAL,INTENT(OUT)             :: WM(1:nVar,1:nGPs)
REAL,INTENT(OUT)             :: WP(1:nVar,1:nGPs)
REAL,INTENT(OUT), OPTIONAL   :: slope_output(1:nVar)
REAL,INTENT(IN), OPTIONAL    :: tau_input
REAL,INTENT(IN), OPTIONAL    :: theta_input
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
!*Default values !*Be careful, once changed they stay changed
REAL, PARAMETER :: theta_default=1.3 
REAL, PARAMETER :: tau_default  =0.5 

REAL             :: theta 
REAL             :: tau   

REAL             :: slope, a, b, c
INTEGER          :: iVar, iGP
REAL             :: w1, w2, w3, rr, phi, rr1, SBM
!-------------------------------------------------------------------------------!

IF (PRESENT(tau_input)) THEN
  tau=tau_input
ELSE
  tau=tau_default
END IF

IF (PRESENT(theta_input)) THEN
  theta=theta_input
ELSE
  theta=theta_default
END IF


DO iGP=1,nGPs
  DO iVar = 1, nVar
    a    = Q(iVar,+1) - Q(iVar,+0)
    b    = Q(iVar,+0) - Q(iVar,-1)

    w1=Q(iVar,-1)
    w2=Q(iVar,+0)
    w3=Q(iVar,+1)



    if(ABS(w2-w1).LT.1e-12)then 
      rr=-1.0
    else
      rr=(w3-w2)/(w2-w1)
    endif

    if(rr .LE. 0.0)then 
      phi=0.0
    elseif(rr .LE. 1.0)then 
      phi=MIN(rr*theta,1.0+tau*(rr-1.0))
    else
      rr1=1.0/rr
      phi=rr*MIN(rr1*theta,1.0+tau*(rr1-1.0))
    endif


    SBM=phi*(w2-w1)


    slope = SBM

    WM(iVar,iGP) = Q(iVar,0) - 0.5*slope
    WP(iVar,iGP) = Q(iVar,0) + 0.5*slope

    IF (PRESENT(slope_output)) THEN
      slope_output(iVar)=slope/dx
    END IF

  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE SBMMUSCL
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE THETAMUSCL(Q,WM,WP,dx,slope_output,theta_input)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)              :: dx
REAL,INTENT(IN)              :: Q(1:nVar,-nGhosts:nGhosts)
REAL,INTENT(OUT)             :: WM(1:nVar,1:nGPs)
REAL,INTENT(OUT)             :: WP(1:nVar,1:nGPs)
REAL,INTENT(OUT), OPTIONAL   :: slope_output(1:nVar)
REAL,INTENT(IN), OPTIONAL    :: theta_input
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
!*Default values !*Be careful, once changed they stay changed
REAL, PARAMETER :: theta_default=1.3 

REAL             :: theta 
REAL             :: tau   

REAL             :: slope, a, b, c
INTEGER          :: iVar, iGP
REAL             :: w1, w2, w3
!-------------------------------------------------------------------------------!

IF (PRESENT(theta_input)) THEN
  theta=theta_input
ELSE
  theta=theta_default
END IF


DO iGP=1,nGPs
  DO iVar = 1, nVar
    a    = Q(iVar,+1) - Q(iVar,+0)
    b    = 0.5*(Q(iVar,+1) - Q(iVar,-1))
    c    = Q(iVar,+0) - Q(iVar,-1)


    slope=MINMOD_3_INPUTS(a*theta,b,c*theta)

    WM(iVar,iGP) = Q(iVar,0) - 0.5*slope
    WP(iVar,iGP) = Q(iVar,0) + 0.5*slope

    IF (PRESENT(slope_output)) THEN
      slope_output(iVar)=slope/dx
    END IF

  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE THETAMUSCL
!===============================================================================!
!
!
!
!===============================================================================!
FUNCTION MINMOD(x,y)
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: x, y
REAL            :: MINMOD
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!

MINMOD = 0.5*(SIGN(1.0,x) + SIGN(1.0,y))*MIN(ABS(x),ABS(y))

!-------------------------------------------------------------------------------!
END FUNCTION MINMOD
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WENO_XDIR(V,WM,WP,dx,WhichReconstruction,RX,LX,RY,LY)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructedVariable
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: dx
REAL,INTENT(IN)    :: V(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts)
REAL,INTENT(OUT)   :: WM(1:nVar,1:nGPs)
REAL,INTENT(OUT)   :: WP(1:nVar,1:nGPs)
INTEGER,INTENT(IN) :: WhichReconstruction
REAL,INTENT(IN)    :: RX(1:nVar,1:nVar)
REAL,INTENT(IN)    :: LX(1:nVar,1:nVar)
REAL,INTENT(IN)    :: RY(1:nVar,1:nVar)
REAL,INTENT(IN)    :: LY(1:nVar,1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: W(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts)
REAL               :: VtempM(1:nVar,-nGhosts:nGhosts)
REAL               :: VtempP(1:nVar,-nGhosts:nGhosts)
INTEGER            :: iVar, ii, jj, kk, ll, iGP
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

!*1) Multiplication by Lx
SELECT CASE(ReconstructedVariable)
  CASE(0,2)
    W(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts)=V(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts)
  CASE(1)
    DO kk=-nGhosts,nGhosts
      DO ll=-nGhosts,nGhosts
        W(1:nVar,kk,ll)=MATMUL(LX(1:nVar,1:nVar),V(1:nVar,kk,ll))
      END DO
    END DO

  CASE DEFAULT
    ErrorMessage = "Error in WENO_XDIR. ReconstructedVariable not available"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!*2) First sweep
DO iVar=1,nVar
  DO jj=-nGhosts,nGhosts
    SELECT CASE(WhichReconstruction)
      CASE(3)
        CALL WENO3_FirstSweep(&
          W(iVar,-nGhosts:nGhosts,jj),VtempM(iVar,jj),VtempP(iVar,jj))
      CASE(4,5)
            CALL WENO5_FirstSweep(&
              W(iVar,-nGhosts:nGhosts,jj),VtempM(iVar,jj),VtempP(iVar,jj))
      CASE(7)
            CALL WENO7_FirstSweep(&
              W(iVar,-nGhosts:nGhosts,jj),VtempM(iVar,jj),VtempP(iVar,jj))
      CASE DEFAULT
        ErrorMessage = "Reconstruction not implemented"
        WRITE(*,*) ErrorMessage
        STOP
    END SELECT
  END DO
END DO

!*3) Multiplication by Rx and, later multiplication by Ly
SELECT CASE(ReconstructedVariable)
  CASE(0,2)
    !*Nothing
  CASE(1)
    DO jj=-nGhosts,nGhosts
      VtempM(1:nVar,jj)=MATMUL(LY,MATMUL(RX,VtempM(1:nVar,jj)))
      VtempP(1:nVar,jj)=MATMUL(LY,MATMUL(RX,VtempP(1:nVar,jj)))
    END DO
  CASE DEFAULT
    ErrorMessage = "Error in WENO_XDIR. ReconstructedVariable not available"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!*4) Second sweep
DO iVar=1,nVar
  SELECT CASE(WhichReconstruction)
    CASE(3)
      CALL WENO3_SecondSweep(VtempM(iVar,-nGhosts:nGhosts),WM(iVar,1:nGPs))
      CALL WENO3_SecondSweep(VtempP(iVar,-nGhosts:nGhosts),WP(iVar,1:nGPs))
    CASE(4,5)
      CALL WENO5_SecondSweep(VtempM(iVar,-nGhosts:nGhosts),WM(iVar,1:nGPs))
      CALL WENO5_SecondSweep(VtempP(iVar,-nGhosts:nGhosts),WP(iVar,1:nGPs))
    CASE(7)
      CALL WENO7_SecondSweep(VtempM(iVar,-nGhosts:nGhosts),WM(iVar,1:nGPs))
      CALL WENO7_SecondSweep(VtempP(iVar,-nGhosts:nGhosts),WP(iVar,1:nGPs))
    CASE DEFAULT
      ErrorMessage = "Reconstruction not implemented in WENO_XDIR"
      WRITE(*,*) ErrorMessage
      STOP
  END SELECT
END DO

!*5) Multiplication by RY
SELECT CASE(ReconstructedVariable)
  CASE(0,2)
    !*Nothing
  CASE(1)
    DO iGP=1,nGPs
      WM(1:nVar,iGP)=MATMUL(RY,WM(1:nVar,iGP))
      WP(1:nVar,iGP)=MATMUL(RY,WP(1:nVar,iGP))
    END DO
  CASE DEFAULT
    ErrorMessage = "Error in WENO_XDIR. ReconstructedVariable not available"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!-------------------------------------------------------------------------------!
END SUBROUTINE WENO_XDIR
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WENO_YDIR(V,WM,WP,dy,WhichReconstruction,RX,LX,RY,LY)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructedVariable
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: dy
REAL,INTENT(IN)    :: V(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts)
REAL,INTENT(OUT)   :: WM(1:nVar,1:nGPs)
REAL,INTENT(OUT)   :: WP(1:nVar,1:nGPs)
INTEGER,INTENT(IN) :: WhichReconstruction
REAL,INTENT(IN)    :: RX(1:nVar,1:nVar)
REAL,INTENT(IN)    :: LX(1:nVar,1:nVar)
REAL,INTENT(IN)    :: RY(1:nVar,1:nVar)
REAL,INTENT(IN)    :: LY(1:nVar,1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: W(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts)
REAL               :: VtempM(1:nVar,-nGhosts:nGhosts)
REAL               :: VtempP(1:nVar,-nGhosts:nGhosts)
INTEGER            :: iVar, ii, jj, kk, ll, iGP
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

!*1) Multiplication by Ly (it was Lx for WENOX)
SELECT CASE(ReconstructedVariable)
  CASE(0,2)
    W(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts)=V(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts)
  CASE(1)
    DO kk=-nGhosts,nGhosts
      DO ll=-nGhosts,nGhosts
        W(1:nVar,kk,ll)=MATMUL(LY(1:nVar,1:nVar),V(1:nVar,kk,ll))
      END DO
    END DO

  CASE DEFAULT
    ErrorMessage = "Error in WENO_YDIR. ReconstructedVariable not available"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!*2) First sweep
DO iVar=1,nVar
  DO ii=-nGhosts,nGhosts
    SELECT CASE(WhichReconstruction)
      CASE(3)
        CALL WENO3_FirstSweep(&
          W(iVar,ii,-nGhosts:nGhosts),VtempM(iVar,ii),VtempP(iVar,ii))
      CASE(4,5)
        CALL WENO5_FirstSweep(&
          W(iVar,ii,-nGhosts:nGhosts),VtempM(iVar,ii),VtempP(iVar,ii))
      CASE(7)
            CALL WENO7_FirstSweep(&
              W(iVar,ii,-nGhosts:nGhosts),VtempM(iVar,ii),VtempP(iVar,ii))
      CASE DEFAULT
        ErrorMessage = "Reconstruction not implemented"
        WRITE(*,*) ErrorMessage
        STOP
    END SELECT
  END DO
END DO


!*3) Multiplication by Ry and, later multiplication by Lx (NB: It was Rx and Ly)
SELECT CASE(ReconstructedVariable)
  CASE(0,2)
    !*Nothing
  CASE(1)
    DO ii=-nGhosts,nGhosts
      VtempM(1:nVar,ii)=MATMUL(LX,MATMUL(RY,VtempM(1:nVar,ii)))
      VtempP(1:nVar,ii)=MATMUL(LX,MATMUL(RY,VtempP(1:nVar,ii)))
    END DO
  CASE DEFAULT
    ErrorMessage = "Error in WENO_YDIR. ReconstructedVariable not available"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!*4) Second sweep
DO iVar=1,nVar
  SELECT CASE(WhichReconstruction)
    CASE(3)
      CALL WENO3_SecondSweep(VtempM(iVar,-nGhosts:nGhosts),WM(iVar,1:nGPs))
      CALL WENO3_SecondSweep(VtempP(iVar,-nGhosts:nGhosts),WP(iVar,1:nGPs))
    CASE(4,5)
      CALL WENO5_SecondSweep(VtempM(iVar,-nGhosts:nGhosts),WM(iVar,1:nGPs))
      CALL WENO5_SecondSweep(VtempP(iVar,-nGhosts:nGhosts),WP(iVar,1:nGPs))
    CASE(7)
      CALL WENO7_SecondSweep(VtempM(iVar,-nGhosts:nGhosts),WM(iVar,1:nGPs))
      CALL WENO7_SecondSweep(VtempP(iVar,-nGhosts:nGhosts),WP(iVar,1:nGPs))
    CASE DEFAULT
      ErrorMessage = "Reconstruction not implemented in WENO_YDIR"
      WRITE(*,*) ErrorMessage
      STOP
  END SELECT
END DO


!*5) Multiplication by RX (it was RY)
SELECT CASE(ReconstructedVariable)
  CASE(0,2)
    !*Nothing
  CASE(1)
    DO iGP=1,nGPs
      WM(1:nVar,iGP)=MATMUL(RX,WM(1:nVar,iGP))
      WP(1:nVar,iGP)=MATMUL(RX,WP(1:nVar,iGP))
    END DO
  CASE DEFAULT
    ErrorMessage = "Error in WENO_YDIR. ReconstructedVariable not available"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE WENO_YDIR
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WENO_In_Quadrature_Points(V,W_in_qp,WhichReconstruction,RX,LX,RY,LY)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructedVariable
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: V(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts)
REAL,INTENT(OUT)   :: W_in_qp(1:nVar,1:nGPs,1:nGPs)
INTEGER,INTENT(IN) :: WhichReconstruction
REAL,INTENT(IN)    :: RX(1:nVar,1:nVar)
REAL,INTENT(IN)    :: LX(1:nVar,1:nVar)
REAL,INTENT(IN)    :: RY(1:nVar,1:nVar)
REAL,INTENT(IN)    :: LY(1:nVar,1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: W(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts)
INTEGER            :: iVar, ii, jj, kk, ll, iGP, jGP, indi, indj
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: Vtemp(1:nVar,-nGhosts:nGhosts,1:nGPs)
REAL               :: Vtemp2(1:nVar,1:nGPs,1:nGPs)
!-------------------------------------------------------------------------------!

!*1) Multiplication by Lx
SELECT CASE(ReconstructedVariable)
  CASE(0,2)
    W(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts)=V(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts)
  CASE(1)
    DO kk=-nGhosts,nGhosts
      DO ll=-nGhosts,nGhosts
        W(1:nVar,kk,ll)=MATMUL(LX(1:nVar,1:nVar),V(1:nVar,kk,ll))
      END DO
    END DO

  CASE DEFAULT
    ErrorMessage = "Error in WENO_In_Quadrature_Points. ReconstructedVariable not available"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!*2) First sweep, X direction
DO iVar=1,nVar
  DO jj=-nGhosts,nGhosts
    SELECT CASE(WhichReconstruction)
      CASE(3)
        CALL WENO3_SecondSweep( W(iVar,-nGhosts:nGhosts,jj), Vtemp(iVar,jj,1:nGPs) )
      CASE(4,5)
        CALL WENO5_SecondSweep( W(iVar,-nGhosts:nGhosts,jj), Vtemp(iVar,jj,1:nGPs) )
      CASE(7)
        CALL WENO7_SecondSweep( W(iVar,-nGhosts:nGhosts,jj), Vtemp(iVar,jj,1:nGPs) )
      CASE DEFAULT
        ErrorMessage = "Reconstruction not implemented"
        WRITE(*,*) ErrorMessage
        STOP
    END SELECT
  END DO
END DO

!*3) Multiplication by Rx and, later multiplication by Ly
SELECT CASE(ReconstructedVariable)
  CASE(0,2)
    !*Nothing
  CASE(1)
    DO jj=-nGhosts,nGhosts
      DO iGP=1,nGPs
        Vtemp(1:nVar,jj,iGP)=MATMUL(LY,MATMUL(RX,Vtemp(1:nVar,jj,iGP)))
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Error in WENO_In_Quadrature_Points. ReconstructedVariable not available"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!*4) Second sweep, Y direction
DO iGP=1,nGPs
  DO iVar=1,nVar
    SELECT CASE(WhichReconstruction)
      CASE(3)
        CALL WENO3_SecondSweep( Vtemp(iVar,-nGhosts:nGhosts,iGP) , Vtemp2(iVar,iGP,1:nGPs) )
      CASE(4,5)
        CALL WENO5_SecondSweep( Vtemp(iVar,-nGhosts:nGhosts,iGP) , Vtemp2(iVar,iGP,1:nGPs) )
      CASE(7)
        CALL WENO7_SecondSweep( Vtemp(iVar,-nGhosts:nGhosts,iGP) , Vtemp2(iVar,iGP,1:nGPs) )
      CASE DEFAULT
        ErrorMessage = "Reconstruction not implemented in WENO_In_Quadrature_Points"
        WRITE(*,*) ErrorMessage
        STOP
    END SELECT
  END DO
END DO

!*5) Multiplication by RY
SELECT CASE(ReconstructedVariable)
  CASE(0,2)
    !*Nothing
  CASE(1)
    DO jGP=1,nGPs
      DO iGP=1,nGPs
        Vtemp2(1:nVar,iGP,jGP)=MATMUL(RY,Vtemp2(1:nVar,iGP,jGP))
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Error in WENO_In_Quadrature_Points. ReconstructedVariable not available"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

W_in_qp(1:nVar,1:nGPs,1:nGPs)=Vtemp2(1:nVar,1:nGPs,1:nGPs)

!-------------------------------------------------------------------------------!
END SUBROUTINE WENO_In_Quadrature_Points
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WENO3_FirstSweep(Q,WM,WP)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
REAL,INTENT(OUT) :: WM
REAL,INTENT(OUT) :: WP
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: alpha1, alpha2
REAL             :: beta1, beta2
REAL             :: gamma1, gamma2
REAL             :: omega1, omega2
REAL             :: W1, W2
!-------------------------------------------------------------------------------!

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
beta1 = (Q(-1) - Q(+0))**2.0
beta2 = (Q(+0) - Q(+1))**2.0


!------------------------------!
! WM: x_{i-1/2}                !
!------------------------------!

! Linear Weights
gamma1 = 2.0/3.0
gamma2 = 1.0/3.0

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2)
omega2 = alpha2/(alpha1 + alpha2)

W1 = 0.5*(    Q(-1) + Q(+0))
W2 = 0.5*(3.0*Q(+0) - Q(+1))
WM = omega1*W1 + omega2*W2


!------------------------------!
! WP: x_{i+1/2}                !
!------------------------------!

! Linear Weights
gamma1 = 1.0/3.0
gamma2 = 2.0/3.0

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2)
omega2 = alpha2/(alpha1 + alpha2)

! Reconstructed Polynomial
W1 = 0.5*(-Q(-1) + 3.0*Q(+0))
W2 = 0.5*( Q(+0) +     Q(+1))
WP = omega1*W1 + omega2*W2

!-------------------------------------------------------------------------------!
END SUBROUTINE WENO3_FirstSweep
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WENO3_SecondSweep(Q,W)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
REAL,INTENT(OUT) :: W(1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: alpha1, alpha2
REAL             :: beta1, beta2
REAL             :: gamma1, gamma2
REAL             :: omega1, omega2
REAL             :: W1, W2
!-------------------------------------------------------------------------------!

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
beta1 = (Q(-1) - Q(+0))**2.0
beta2 = (Q(+0) - Q(+1))**2.0


!------------------------------!
! Point: x_{j-1/(2*sqrt(3))}   !
!------------------------------!

! Linear Weights
gamma1 = 0.5
gamma2 = 0.5

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2)
omega2 = alpha2/(alpha1 + alpha2)

! Reconstructed Polynomial
W1   = (1.0/6.0)*(SQRT(3.0)*Q(-1) + 6.0*Q(+0) - SQRT(3.0)*Q(+0))
W2   = (1.0/6.0)*(SQRT(3.0)*Q(+0) + 6.0*Q(+0) - SQRT(3.0)*Q(+1))
W(1) = omega1*W1 + omega2*W2


!------------------------------!
! Point: x_{j+1/(2*sqrt(3))}   !
!------------------------------!

! Linear Weights
gamma1 = 0.5
gamma2 = 0.5

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2)
omega2 = alpha2/(alpha1 + alpha2)

! Reconstructed Polynomial
W1   = (1.0/6.0)*(-SQRT(3.0)*Q(-1) + 6.0*Q(+0) + SQRT(3.0)*Q(+0))
W2   = (1.0/6.0)*(-SQRT(3.0)*Q(+0) + 6.0*Q(+0) + SQRT(3.0)*Q(+1))
W(2) = omega1*W1 + omega2*W2

!-------------------------------------------------------------------------------!
END SUBROUTINE WENO3_SecondSweep
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WENO5_FirstSweep(Q,WM,WP)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
REAL,INTENT(OUT) :: WM
REAL,INTENT(OUT) :: WP
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: alpha1, alpha2, alpha3
REAL             :: beta1,  beta2,  beta3
REAL             :: gamma1, gamma2, gamma3
REAL             :: omega1, omega2, omega3
REAL             :: W1, W2, W3
!-------------------------------------------------------------------------------!

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
beta1 = (1.0/3.0)*( 4.0*Q(-2)*Q(-2) - 19.0*Q(-2)*Q(-1) + 25.0*Q(-1)*Q(-1) &
                 + 11.0*Q(-2)*Q(+0) - 31.0*Q(-1)*Q(+0) + 10.0*Q(+0)*Q(+0))
beta2 = (1.0/3.0)*( 4.0*Q(-1)*Q(-1) - 13.0*Q(-1)*Q(+0) + 13.0*Q(+0)*Q(+0) &
                 +  5.0*Q(-1)*Q(+1) - 13.0*Q(+0)*Q(+1) +  4.0*Q(+1)*Q(+1))
beta3 = (1.0/3.0)*(10.0*Q(+0)*Q(+0) - 31.0*Q(+0)*Q(+1) + 25.0*Q(+1)*Q(+1) &
                 + 11.0*Q(+0)*Q(+2) - 19.0*Q(+1)*Q(+2) +  4.0*Q(+2)*Q(+2))


!------------------------------!
! WM: x_{i-1/2}                !
!------------------------------!

! Linear Weights
gamma1 = 3.0/10.0
gamma2 = 3.0/5.0
gamma3 = 1.0/10.0

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

W1 = (1.0/6.0)*(    -Q(-2) + 5.0*Q(-1) + 2.0*Q(+0))
W2 = (1.0/6.0)*( 2.0*Q(-1) + 5.0*Q(+0) -     Q(+1))
W3 = (1.0/6.0)*(11.0*Q(+0) - 7.0*Q(+1) + 2.0*Q(+2))
WM = omega1*W1 + omega2*W2 + omega3*W3


!------------------------------!
! WP: x_{i+1/2}                !
!------------------------------!

! Linear Weights
gamma1 = 1.0/10.0
gamma2 = 3.0/5.0
gamma3 = 3.0/10.0

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

W1 = (1.0/6.0)*(2.0*Q(-2) - 7.0*Q(-1) + 11.0*Q(+0))
W2 = (1.0/6.0)*(   -Q(-1) + 5.0*Q(+0) +  2.0*Q(+1))
W3 = (1.0/6.0)*(2.0*Q(+0) + 5.0*Q(+1) -      Q(+2))
WP = omega1*W1 + omega2*W2 + omega3*W3

!-------------------------------------------------------------------------------!
END SUBROUTINE WENO5_FirstSweep
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WENO5_SecondSweep2nGPs(Q,W)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
REAL,INTENT(OUT) :: W(1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: alpha1, alpha2, alpha3
REAL             :: beta1,  beta2,  beta3
REAL             :: gamma1, gamma2, gamma3
REAL             :: omega1, omega2, omega3
REAL             :: W1, W2, W3
!-------------------------------------------------------------------------------!

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
beta1 = (1.0/3.0)*( 4.0*Q(-2)*Q(-2) - 19.0*Q(-2)*Q(-1) + 25.0*Q(-1)*Q(-1) &
                 + 11.0*Q(-2)*Q(+0) - 31.0*Q(-1)*Q(+0) + 10.0*Q(+0)*Q(+0))
beta2 = (1.0/3.0)*( 4.0*Q(-1)*Q(-1) - 13.0*Q(-1)*Q(+0) + 13.0*Q(+0)*Q(+0) &
                 +  5.0*Q(-1)*Q(+1) - 13.0*Q(+0)*Q(+1) +  4.0*Q(+1)*Q(+1))
beta3 = (1.0/3.0)*(10.0*Q(+0)*Q(+0) - 31.0*Q(+0)*Q(+1) + 25.0*Q(+1)*Q(+1) &
                 + 11.0*Q(+0)*Q(+2) - 19.0*Q(+1)*Q(+2) +  4.0*Q(+2)*Q(+2))


!------------------------------!
! Point: x_{j-1/(2*sqrt(3))}   !
!------------------------------!

! Linear Weights
gamma1 = (210.0 + SQRT(3.0))/1080.0
gamma2 = 11.0/18.0
gamma3 = (210.0 - SQRT(3.0))/1080.0

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial
W1     = -     SQRT(3.0)*Q(-2) &
         + 4.0*SQRT(3.0)*Q(-1) &
         +          12.0*Q(+0) &
         - 3.0*SQRT(3.0)*Q(+0)
W2     = + 1.0*SQRT(3.0)*Q(-1) &
         +          12.0*Q(+0) &
         - 1.0*SQRT(3.0)*Q(+1)
W3     = +          12.0*Q(+0) &
         + 3.0*SQRT(3.0)*Q(+0) &
         - 4.0*SQRT(3.0)*Q(+1) &
         +     SQRT(3.0)*Q(+2)
W(1)   = (omega1*W1 + omega2*W2 + omega3*W3)/12.0


!------------------------------!
! Point: x_{j+1/(2*sqrt(3))}   !
!------------------------------!

! Linear Weights
gamma1 = (210.0 - SQRT(3.0))/1080.0
gamma2 = 11.0/18.0
gamma3 = (210.0 + SQRT(3.0))/1080.0

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial
W1     = + 1.0*SQRT(3.0)*Q(-2) &
         - 4.0*SQRT(3.0)*Q(-1) &
         +          12.0*Q(+0) &
         + 3.0*SQRT(3.0)*Q(+0)
W2     = - 1.0*SQRT(3.0)*Q(-1) &
         +          12.0*Q(+0) &
         +     SQRT(3.0)*Q(+1)
W3     = +          12.0*Q(+0) &
         - 3.0*SQRT(3.0)*Q(+0) &
         + 4.0*SQRT(3.0)*Q(+1) &
         -     SQRT(3.0)*Q(+2)
W(2)   = (omega1*W1 + omega2*W2 + omega3*W3)/12.0

!-------------------------------------------------------------------------------!
END SUBROUTINE WENO5_SecondSweep2nGPs
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WENO5_SecondSweep3nGPs(Q,W)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
REAL,INTENT(OUT) :: W(1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: alpha1, alpha2, alpha3
REAL             :: beta1,  beta2,  beta3
REAL             :: gamma1, gamma2, gamma3
REAL             :: omega1, omega2, omega3
REAL             :: W1, W2, W3
!-------------------------------------------------------------------------------!

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
beta1 = (1.0/3.0)*( 4.0*Q(-2)*Q(-2) - 19.0*Q(-2)*Q(-1) + 25.0*Q(-1)*Q(-1) &
                 + 11.0*Q(-2)*Q(+0) - 31.0*Q(-1)*Q(+0) + 10.0*Q(+0)*Q(+0))
beta2 = (1.0/3.0)*( 4.0*Q(-1)*Q(-1) - 13.0*Q(-1)*Q(+0) + 13.0*Q(+0)*Q(+0) &
                 +  5.0*Q(-1)*Q(+1) - 13.0*Q(+0)*Q(+1) +  4.0*Q(+1)*Q(+1))
beta3 = (1.0/3.0)*(10.0*Q(+0)*Q(+0) - 31.0*Q(+0)*Q(+1) + 25.0*Q(+1)*Q(+1) &
                 + 11.0*Q(+0)*Q(+2) - 19.0*Q(+1)*Q(+2) +  4.0*Q(+2)*Q(+2))

!------------------------------!
! Point: x_{j-1/2*sqrt(3/5)}   !
!------------------------------!

! Linear Weights
gamma1 = (71.0*SQRT(15.0))/5240.0 + 126.0/655.0
gamma2 = 403.0/655.0
gamma3 = -(71.0*SQRT(15.0))/5240.0 + 126.0/655.0

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial
W1     = (1.0/30.0 - SQRT(15.)/20.0)*Q(-2) &
         + (SQRT(15.)/5. - 1./15.)*Q(-1) &
         + (31./30. - (3.*SQRT(15.))/20.)*Q(+0)
W2     = + (SQRT(15.)/20. + 1./30.)*Q(-1) &
         +         14./15.*Q(+0) &
         + (1./30. - SQRT(15.)/20.)*Q(+1)
W3     = + (3.*SQRT(15.)/20. + 31./30.)*Q(+0) &
         - (SQRT(15.)/5. + 1./15.)*Q(+1) &
         + (SQRT(15.)/20. + 1./30.)*Q(+2)
W(1)   = (omega1*W1 + omega2*W2 + omega3*W3)


!------------------------------!
! Point: x_{j}                 !
!------------------------------!

! Linear Weights
gamma1 = -9./80.
gamma2 = 49./40.
gamma3 = -9./80.

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial
W1     = -1./24. *Q(-2) &
         +1./12. *Q(-1) &
         +23./24. *Q(+0) 
W2     = - 1./24. *Q(-1) &
         + 13./12.*Q(+0) &
         - 1./24. *Q(+1)
W3     = + 23./24.*Q(+0) &
         + 1. /12.*Q(+1) &
         - 1./24. *Q(+2)

W(2)   = (omega1*W1 + omega2*W2 + omega3*W3)



!------------------------------!
! Point: x_{j+1/2*sqrt(3/5)}   !
!------------------------------!

! Linear Weights
gamma1 = -(71.0*SQRT(15.0))/5240.0 + 126.0/655.0
gamma2 = 403.0/655.0
gamma3 = +(71.0*SQRT(15.0))/5240.0 + 126.0/655.0

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial
W1     = (1.0/30.0 + SQRT(15.)/20.0)*Q(-2) &
         - (SQRT(15.)/5. + 1./15.)*Q(-1) &
         + (31./30. + (3.*SQRT(15.))/20.)*Q(+0)
W2     = + (-SQRT(15.)/20. + 1./30.)*Q(-1) &
         +         14./15.*Q(+0) &
         + (1./30. + SQRT(15.)/20.)*Q(+1)
W3     = + (-3.*SQRT(15.)/20. + 31./30.)*Q(+0) &
         + (SQRT(15.)/5. - 1./15.)*Q(+1) &
         + (-SQRT(15.)/20. + 1./30.)*Q(+2)
W(3)   = (omega1*W1 + omega2*W2 + omega3*W3)

!-------------------------------------------------------------------------------!
END SUBROUTINE WENO5_SecondSweep3nGPs
!===============================================================================!
!
!!
!
!===============================================================================!
SUBROUTINE WENO5_SecondSweep(Q,W) !4nGPS
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
REAL,INTENT(OUT) :: W(1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: alpha1, alpha2, alpha3
REAL             :: beta1,  beta2,  beta3
REAL             :: gamma1, gamma2, gamma3
REAL             :: omega1, omega2, omega3
REAL             :: W1, W2, W3
!-------------------------------------------------------------------------------!

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
beta1 = (1.0/3.0)*( 4.0*Q(-2)*Q(-2) - 19.0*Q(-2)*Q(-1) + 25.0*Q(-1)*Q(-1) &
                 + 11.0*Q(-2)*Q(+0) - 31.0*Q(-1)*Q(+0) + 10.0*Q(+0)*Q(+0))
beta2 = (1.0/3.0)*( 4.0*Q(-1)*Q(-1) - 13.0*Q(-1)*Q(+0) + 13.0*Q(+0)*Q(+0) &
                 +  5.0*Q(-1)*Q(+1) - 13.0*Q(+0)*Q(+1) +  4.0*Q(+1)*Q(+1))
beta3 = (1.0/3.0)*(10.0*Q(+0)*Q(+0) - 31.0*Q(+0)*Q(+1) + 25.0*Q(+1)*Q(+1) &
                 + 11.0*Q(+0)*Q(+2) - 19.0*Q(+1)*Q(+2) +  4.0*Q(+2)*Q(+2))

!--------------------------------------------!
! Point: x_{j-1/2*sqrt(3/7+2/7*sqrt(6/5))}   !
!--------------------------------------------!

! Linear Weights
gamma1 =0.2658420974778319
gamma2 =0.6112504900322486
gamma3 =0.1229074124899195

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial
W1 = -0.1642562761719537 * Q(-2)+0.7590807081409336 * Q(-1)+0.4051755680310201 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)
W2 = +0.0000000000000000 * Q(-2)+0.2663118796250726 * Q(-1)+0.8979443965468811 * Q(0)-0.1642562761719537 * Q(1)+0.0000000000000000 * Q(2)
W3 = +0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+1.6968800354220990 * Q(0)-0.9631919150471715 * Q(1)+0.2663118796250726 * Q(2)
W(1)   = (omega1*W1 + omega2*W2 + omega3*W3)

!--------------------------------------------!
! Point: x_{j-1/2*sqrt(3/7-2/7*sqrt(6/5))}   !
!--------------------------------------------!

! Linear Weights

gamma1 =0.1281641584355059
gamma2 =0.5219691498503563
gamma3 =0.3498666917141379

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)


W1 = -0.1122135388132497 * Q(-2)+0.3944175994189276 * Q(-1)+0.7177959393943222 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)
W2 = +0.0000000000000000 * Q(-2)+0.0577769829791784 * Q(-1)+1.0544365558340714 * Q(0)-0.1122135388132497 * Q(1)+0.0000000000000000 * Q(2)
W3 = +0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+1.2277675047716066 * Q(0)-0.2855444877507849 * Q(1)+0.0577769829791784 * Q(2)

W(2)   = (omega1*W1 + omega2*W2 + omega3*W3)


!--------------------------------------------!
! Point: x_{j+1/2*sqrt(3/7-2/7*sqrt(6/5))}   !
!--------------------------------------------!

! Linear Weights
gamma1 =0.3498666917141379
gamma2 =0.5219691498503563
gamma3 =0.1281641584355059

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial

W1 = +0.0577769829791784 * Q(-2)-0.2855444877507849 * Q(-1)+1.2277675047716066 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)
W2 = +0.0000000000000000 * Q(-2)-0.1122135388132497 * Q(-1)+1.0544365558340714 * Q(0)+0.0577769829791784 * Q(1)+0.0000000000000000 * Q(2)
W3 = +0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+0.7177959393943222 * Q(0)+0.3944175994189276 * Q(1)-0.1122135388132497 * Q(2)

W(3)   = (omega1*W1 + omega2*W2 + omega3*W3)



!--------------------------------------------!
! Point: x_{j+1/2*sqrt(3/7+2/7*sqrt(6/5))}   !
!--------------------------------------------!

! Linear Weights
gamma1 =0.1229074124899195
gamma2 =0.6112504900322486
gamma3 =0.2658420974778319


alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial
W1 = +0.2663118796250726 * Q(-2)-0.9631919150471715 * Q(-1)+1.6968800354220990 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)
W2 = +0.0000000000000000 * Q(-2)-0.1642562761719537 * Q(-1)+0.8979443965468811 * Q(0)+0.2663118796250726 * Q(1)+0.0000000000000000 * Q(2)
W3 = +0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+0.4051755680310201 * Q(0)+0.7590807081409336 * Q(1)-0.1642562761719537 * Q(2)
W(4)   = (omega1*W1 + omega2*W2 + omega3*W3)

!-------------------------------------------------------------------------------!
END SUBROUTINE WENO5_SecondSweep
!===============================================================================!
!
!
!
!
!===============================================================================!
SUBROUTINE WENO7_FirstSweep(Q,WM,WP)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
REAL,INTENT(OUT) :: WM
REAL,INTENT(OUT) :: WP
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: alpha1, alpha2, alpha3, alpha4
REAL             :: beta1,  beta2,  beta3,  beta4
REAL             :: gamma1, gamma2, gamma3, gamma4
REAL             :: omega1, omega2, omega3, omega4
REAL             :: W1, W2, W3, W4
!-------------------------------------------------------------------------------!

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
! beta1 = +(547./240.*Q(-3)*Q(-3))+(-647./80.*Q(-3)*Q(-2))+(2321./240.*Q(-3)*Q(-1))+(-309./80.*Q(-3)*Q(0))+(-647./80.*Q(-2)*Q(-3))+(7043./240.*Q(-2)*Q(-2))+(-8623./240.*Q(-2)*Q(-1))+(3521./240.*Q(-2)*Q(0))+(2321./240.*Q(-1)*Q(-3))+(-8623./240.*Q(-1)*Q(-2))+(11003/240*Q(-1)*Q(-1))+(-1567/80*Q(-1)*Q(0))+(-309/80*Q(0)*Q(-3))+(3521/240*Q(0)*Q(-2))+(-1567/80*Q(0)*Q(-1))+(2107/240*Q(0)*Q(0))
! beta2 = +(89/80*Q(-2)*Q(-2))+(-821/240*Q(-2)*Q(-1))+(267/80*Q(-2)*Q(0))+(-247/240*Q(-2)*Q(1))+(-821/240*Q(-1)*Q(-2))+(2843/240*Q(-1)*Q(-1))+(-2983/240*Q(-1)*Q(0))+(961/240*Q(-1)*Q(1))+(267/80*Q(0)*Q(-2))+(-2983/240*Q(0)*Q(-1))+(3443/240*Q(0)*Q(0))+(-1261/240*Q(0)*Q(1))+(-247/240*Q(1)*Q(-2))+(961/240*Q(1)*Q(-1))+(-1261/240*Q(1)*Q(0))+(547/240*Q(1)*Q(1))
! beta3 = +(547/240*Q(-1)*Q(-1))+(-1261/240*Q(-1)*Q(0))+(961/240*Q(-1)*Q(1))+(-247/240*Q(-1)*Q(2))+(-1261/240*Q(0)*Q(-1))+(3443/240*Q(0)*Q(0))+(-2983/240*Q(0)*Q(1))+(267/80*Q(0)*Q(2))+(961/240*Q(1)*Q(-1))+(-2983/240*Q(1)*Q(0))+(2843/240*Q(1)*Q(1))+(-821/240*Q(1)*Q(2))+(-247/240*Q(2)*Q(-1))+(267/80*Q(2)*Q(0))+(-821/240*Q(2)*Q(1))+(89/80*Q(2)*Q(2))
! beta4 = +(2107/240*Q(0)*Q(0))+(-1567/80*Q(0)*Q(1))+(3521/240*Q(0)*Q(2))+(-309/80*Q(0)*Q(3))+(-1567/80*Q(1)*Q(0))+(11003/240*Q(1)*Q(1))+(-8623/240*Q(1)*Q(2))+(2321/240*Q(1)*Q(3))+(3521/240*Q(2)*Q(0))+(-8623/240*Q(2)*Q(1))+(7043/240*Q(2)*Q(2))+(-647/80*Q(2)*Q(3))+(-309/80*Q(3)*Q(0))+(2321/240*Q(3)*Q(1))+(-647/80*Q(3)*Q(2))+(547/240*Q(3)*Q(3))

beta1 = +(2.27916666666666679*Q(-3)*Q(-3))+(-8.08750000000000036*Q(-3)*Q(-2))+(9.67083333333333250*Q(-3)*Q(-1))+(-3.8625000000000000*Q(-3)*Q(0))+(-8.08750000000000036*Q(-2)*Q(-3))+(29.34583333333333499*Q(-2)*Q(-2))+(-35.92916666666666714*Q(-2)*Q(-1))+(14.67083333333333250*Q(-2)*Q(0))+(9.67083333333333250*Q(-1)*Q(-3))+(-35.92916666666666714*Q(-1)*Q(-2))+(45.84583333333333144*Q(-1)*Q(-1))+(-19.5875000000000000*Q(-1)*Q(0))+(-3.8625000000000000*Q(0)*Q(-3))+(14.67083333333333250*Q(0)*Q(-2))+(-19.58749999999999858*Q(0)*Q(-1))+(8.77916666666666679*Q(0)*Q(0))
beta2 = +(1.11250000000000004*Q(-2)*Q(-2))+(-3.42083333333333339*Q(-2)*Q(-1))+(3.33749999999999991*Q(-2)*Q(0))+(-1.02916666666666656*Q(-2)*Q(1))+(-3.42083333333333339*Q(-1)*Q(-2))+(11.84583333333333321*Q(-1)*Q(-1))+(-12.42916666666666714*Q(-1)*Q(0))+(4.00416666666666643*Q(-1)*Q(1))+(3.3375000000000000*Q(0)*Q(-2))+(-12.42916666666666714*Q(0)*Q(-1))+(14.34583333333333321*Q(0)*Q(0))+(-5.25416666666666643*Q(0)*Q(1))+(-1.02916666666666656*Q(1)*Q(-2))+(4.00416666666666643*Q(1)*Q(-1))+(-5.25416666666666643*Q(1)*Q(0))+(2.27916666666666679*Q(1)*Q(1))
beta3 = +(2.27916666666666679*Q(-1)*Q(-1))+(-5.25416666666666643*Q(-1)*Q(0))+(4.00416666666666643*Q(-1)*Q(1))+(-1.02916666666666656*Q(-1)*Q(2))+(-5.25416666666666643*Q(0)*Q(-1))+(14.34583333333333321*Q(0)*Q(0))+(-12.42916666666666714*Q(0)*Q(1))+(3.33749999999999991*Q(0)*Q(2))+(4.00416666666666643*Q(1)*Q(-1))+(-12.42916666666666714*Q(1)*Q(0))+(11.84583333333333321*Q(1)*Q(1))+(-3.42083333333333339*Q(1)*Q(2))+(-1.02916666666666656*Q(2)*Q(-1))+(3.3375000000000000*Q(2)*Q(0))+(-3.42083333333333339*Q(2)*Q(1))+(1.11250000000000004*Q(2)*Q(2))
beta4 = +(8.77916666666666679*Q(0)*Q(0))+(-19.5875000000000000*Q(0)*Q(1))+(14.67083333333333250*Q(0)*Q(2))+(-3.8625000000000000*Q(0)*Q(3))+(-19.5875000000000000*Q(1)*Q(0))+(45.84583333333333144*Q(1)*Q(1))+(-35.92916666666666714*Q(1)*Q(2))+(9.67083333333333250*Q(1)*Q(3))+(14.67083333333333250*Q(2)*Q(0))+(-35.92916666666666714*Q(2)*Q(1))+(29.34583333333333499*Q(2)*Q(2))+(-8.08750000000000036*Q(2)*Q(3))+(-3.86249999999999982*Q(3)*Q(0))+(9.67083333333333250*Q(3)*Q(1))+(-8.08750000000000036*Q(3)*Q(2))+(2.27916666666666679*Q(3)*Q(3))


!------------------------------!
! WM: x_{i-1/2}                !
!------------------------------!

! Linear Weights
gamma1 =0.1142857142857143
gamma2 =0.5142857142857142
gamma3 =0.3428571428571429
gamma4 =0.0285714285714286

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP
alpha4 = gamma4/(WENOEPS + beta4)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3 + alpha4)
omega2 = alpha2/(alpha1 + alpha2 + alpha3 + alpha4)
omega3 = alpha3/(alpha1 + alpha2 + alpha3 + alpha4)
omega4 = alpha4/(alpha1 + alpha2 + alpha3 + alpha4)

W1 = +0.0833333333333333 * Q(-3)-0.4166666666666667 * Q(-2)+1.0833333333333333 * Q(-1)+0.2500000000000000 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)
W2 = +0.0000000000000000 * Q(-3)-0.0833333333333333 * Q(-2)+0.5833333333333334 * Q(-1)+0.5833333333333334 * Q(0)-0.0833333333333333 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)
W3 = +0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.2500000000000000 * Q(-1)+1.0833333333333333 * Q(0)-0.4166666666666667 * Q(1)+0.0833333333333333 * Q(2)+0.0000000000000000 * Q(3)
W4 = +0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+2.0833333333333335 * Q(0)-1.9166666666666667 * Q(1)+1.0833333333333333 * Q(2)-0.2500000000000000 * Q(3)



  WM = omega1*W1 + omega2*W2 + omega3*W3 + omega4*W4


!------------------------------!
! WP: x_{i+1/2}                !
!------------------------------!

! Linear Weights
gamma1 =0.0285714285714286
gamma2 =0.3428571428571429
gamma3 =0.5142857142857142
gamma4 =0.1142857142857143

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP
alpha4 = gamma4/(WENOEPS + beta4)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3 + alpha4)
omega2 = alpha2/(alpha1 + alpha2 + alpha3 + alpha4)
omega3 = alpha3/(alpha1 + alpha2 + alpha3 + alpha4)
omega4 = alpha4/(alpha1 + alpha2 + alpha3 + alpha4)

W1 = -0.2500000000000000 * Q(-3)+1.0833333333333333 * Q(-2)-1.9166666666666667 * Q(-1)+2.0833333333333335 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)
W2 = +0.0000000000000000 * Q(-3)+0.0833333333333333 * Q(-2)-0.4166666666666667 * Q(-1)+1.0833333333333333 * Q(0)+0.2500000000000000 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)
W3 = +0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)-0.0833333333333333 * Q(-1)+0.5833333333333334 * Q(0)+0.5833333333333334 * Q(1)-0.0833333333333333 * Q(2)+0.0000000000000000 * Q(3)
W4 = +0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+0.2500000000000000 * Q(0)+1.0833333333333333 * Q(1)-0.4166666666666667 * Q(2)+0.0833333333333333 * Q(3)


  WP = omega1*W1 + omega2*W2 + omega3*W3 + omega4*W4

!-------------------------------------------------------------------------------!
END SUBROUTINE WENO7_FirstSweep
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WENO7_SecondSweep(Q,W) !4nGPS
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
REAL,INTENT(OUT) :: W(1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: alpha1, alpha2, alpha3, alpha4
REAL             :: beta1,  beta2,  beta3,  beta4
REAL             :: gamma1, gamma2, gamma3, gamma4
REAL             :: omega1, omega2, omega3, omega4
REAL             :: W1, W2, W3, W4
!-------------------------------------------------------------------------------!

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
beta1 = +(2.27916666666666679*Q(-3)*Q(-3))+(-8.08750000000000036*Q(-3)*Q(-2))+(9.67083333333333250*Q(-3)*Q(-1))+(-3.86249999999999982*Q(-3)*Q(0))+(-8.08750000000000036*Q(-2)*Q(-3))+(29.34583333333333499*Q(-2)*Q(-2))+(-35.92916666666666714*Q(-2)*Q(-1))+(14.67083333333333250*Q(-2)*Q(0))+(9.67083333333333250*Q(-1)*Q(-3))+(-35.92916666666666714*Q(-1)*Q(-2))+(45.84583333333333144*Q(-1)*Q(-1))+(-19.58749999999999858*Q(-1)*Q(0))+(-3.86249999999999982*Q(0)*Q(-3))+(14.67083333333333250*Q(0)*Q(-2))+(-19.58749999999999858*Q(0)*Q(-1))+(8.77916666666666679*Q(0)*Q(0))
beta2 = +(1.11250000000000004*Q(-2)*Q(-2))+(-3.42083333333333339*Q(-2)*Q(-1))+(3.33749999999999991*Q(-2)*Q(0))+(-1.02916666666666656*Q(-2)*Q(1))+(-3.42083333333333339*Q(-1)*Q(-2))+(11.84583333333333321*Q(-1)*Q(-1))+(-12.42916666666666714*Q(-1)*Q(0))+(4.00416666666666643*Q(-1)*Q(1))+(3.33749999999999991*Q(0)*Q(-2))+(-12.42916666666666714*Q(0)*Q(-1))+(14.34583333333333321*Q(0)*Q(0))+(-5.25416666666666643*Q(0)*Q(1))+(-1.02916666666666656*Q(1)*Q(-2))+(4.00416666666666643*Q(1)*Q(-1))+(-5.25416666666666643*Q(1)*Q(0))+(2.27916666666666679*Q(1)*Q(1))
beta3 = +(2.27916666666666679*Q(-1)*Q(-1))+(-5.25416666666666643*Q(-1)*Q(0))+(4.00416666666666643*Q(-1)*Q(1))+(-1.02916666666666656*Q(-1)*Q(2))+(-5.25416666666666643*Q(0)*Q(-1))+(14.34583333333333321*Q(0)*Q(0))+(-12.42916666666666714*Q(0)*Q(1))+(3.33749999999999991*Q(0)*Q(2))+(4.00416666666666643*Q(1)*Q(-1))+(-12.42916666666666714*Q(1)*Q(0))+(11.84583333333333321*Q(1)*Q(1))+(-3.42083333333333339*Q(1)*Q(2))+(-1.02916666666666656*Q(2)*Q(-1))+(3.33749999999999991*Q(2)*Q(0))+(-3.42083333333333339*Q(2)*Q(1))+(1.11250000000000004*Q(2)*Q(2))
beta4 = +(8.77916666666666679*Q(0)*Q(0))+(-19.58749999999999858*Q(0)*Q(1))+(14.67083333333333250*Q(0)*Q(2))+(-3.86249999999999982*Q(0)*Q(3))+(-19.58749999999999858*Q(1)*Q(0))+(45.84583333333333144*Q(1)*Q(1))+(-35.92916666666666714*Q(1)*Q(2))+(9.67083333333333250*Q(1)*Q(3))+(14.67083333333333250*Q(2)*Q(0))+(-35.92916666666666714*Q(2)*Q(1))+(29.34583333333333499*Q(2)*Q(2))+(-8.08750000000000036*Q(2)*Q(3))+(-3.86249999999999982*Q(3)*Q(0))+(9.67083333333333250*Q(3)*Q(1))+(-8.08750000000000036*Q(3)*Q(2))+(2.27916666666666679*Q(3)*Q(3))

!--------------------------------------------!
! Point: x_{j-1/2*sqrt(3/7+2/7*sqrt(6/5))}   !
!--------------------------------------------!

! Linear Weights
gamma1 =0.0978973393748262
gamma2 =0.4945893001737974
gamma3 =0.3706410089227926
gamma4 =0.0368723515285839

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP
alpha4 = gamma4/(WENOEPS + beta4)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3 + alpha4)
omega2 = alpha2/(alpha1 + alpha2 + alpha3 + alpha4)
omega3 = alpha3/(alpha1 + alpha2 + alpha3 + alpha4)
omega4 = alpha4/(alpha1 + alpha2 + alpha3 + alpha4)

! Reconstructed Polynomial
W1 = +0.0878583391504589 * Q(-3)-0.4278312936233303 * Q(-2)+1.0226557255923103 * Q(-1)+0.3173172288805612 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)
W2 = +0.0000000000000000 * Q(-3)-0.0763979370214948 * Q(-2)+0.4955056906895569 * Q(-1)+0.6687505854823967 * Q(0)-0.0878583391504589 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)
W3 = +0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.1899139426035779 * Q(-1)+1.1271382076113654 * Q(0)-0.3934500872364380 * Q(1)+0.0763979370214948 * Q(2)+0.0000000000000000 * Q(3)
W4 = +0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+1.8867939780256768 * Q(0)-1.5329337428579051 * Q(1)+0.8360537074358062 * Q(2)-0.1899139426035779 * Q(3)

  W(1)   = (omega1*W1 + omega2*W2 + omega3*W3 +omega4*W4)


!--------------------------------------------!
! Point: x_{j-1/2*sqrt(3/7-2/7*sqrt(6/5))}   !
!--------------------------------------------!

! Linear Weights

gamma1 =0.0422160570229445
gamma2 =0.3488436391325038
gamma3 =0.4308046279496520
gamma4 =0.1781356758948997

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP
alpha4 = gamma4/(WENOEPS + beta4)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3 + alpha4)
omega2 = alpha2/(alpha1 + alpha2 + alpha3 + alpha4)
omega3 = alpha3/(alpha1 + alpha2 + alpha3 + alpha4)
omega4 = alpha4/(alpha1 + alpha2 + alpha3 + alpha4)


W1 = +0.0776175431540304 * Q(-3)-0.3450661682753410 * Q(-2)+0.6272702288810189 * Q(-1)+0.6401783962402917 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)
W2 = +0.0000000000000000 * Q(-3)-0.0345959956592193 * Q(-2)+0.1615649699568364 * Q(-1)+0.9506485688564134 * Q(0)-0.0776175431540304 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)
W3 = +0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.0231809873199591 * Q(-1)+1.1582245428117293 * Q(0)-0.2160015257909077 * Q(1)+0.0345959956592193 * Q(2)+0.0000000000000000 * Q(3)
W4 = +0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+1.2509484920915657 * Q(0)-0.3550874497106621 * Q(1)+0.1273199449390556 * Q(2)-0.0231809873199591 * Q(3)

  W(2)   = (omega1*W1 + omega2*W2 + omega3*W3 +omega4*W4)



!--------------------------------------------!
! Point: x_{j+1/2*sqrt(3/7-2/7*sqrt(6/5))}   !
!--------------------------------------------!

! Linear Weights
gamma1 =0.1781356758948997
gamma2 =0.4308046279496520
gamma3 =0.3488436391325038
gamma4 =0.0422160570229445

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP
alpha4 = gamma4/(WENOEPS + beta4)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3 + alpha4)
omega2 = alpha2/(alpha1 + alpha2 + alpha3 + alpha4)
omega3 = alpha3/(alpha1 + alpha2 + alpha3 + alpha4)
omega4 = alpha4/(alpha1 + alpha2 + alpha3 + alpha4)


! Reconstructed Polynomial

W1 = -0.0231809873199591 * Q(-3)+0.1273199449390556 * Q(-2)-0.3550874497106621 * Q(-1)+1.2509484920915657 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)
W2 = +0.0000000000000000 * Q(-3)+0.0345959956592193 * Q(-2)-0.2160015257909077 * Q(-1)+1.1582245428117293 * Q(0)+0.0231809873199591 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)
W3 = +0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)-0.0776175431540304 * Q(-1)+0.9506485688564134 * Q(0)+0.1615649699568364 * Q(1)-0.0345959956592193 * Q(2)+0.0000000000000000 * Q(3)
W4 = +0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+0.6401783962402917 * Q(0)+0.6272702288810189 * Q(1)-0.3450661682753410 * Q(2)+0.0776175431540304 * Q(3)

  W(3)   = (omega1*W1 + omega2*W2 + omega3*W3 +omega4*W4)



!--------------------------------------------!
! Point: x_{j+1/2*sqrt(3/7+2/7*sqrt(6/5))}   !
!--------------------------------------------!

! Linear Weights
gamma1 =0.0368723515285839
gamma2 =0.3706410089227926
gamma3 =0.4945893001737974
gamma4 =0.0978973393748262

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP
alpha4 = gamma4/(WENOEPS + beta4)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3 + alpha4)
omega2 = alpha2/(alpha1 + alpha2 + alpha3 + alpha4)
omega3 = alpha3/(alpha1 + alpha2 + alpha3 + alpha4)
omega4 = alpha4/(alpha1 + alpha2 + alpha3 + alpha4)

! Reconstructed Polynomial
W1 = -0.1899139426035779 * Q(-3)+0.8360537074358062 * Q(-2)-1.5329337428579051 * Q(-1)+1.8867939780256768 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)
W2 = +0.0000000000000000 * Q(-3)+0.0763979370214948 * Q(-2)-0.3934500872364380 * Q(-1)+1.1271382076113654 * Q(0)+0.1899139426035779 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)
W3 = +0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)-0.0878583391504589 * Q(-1)+0.6687505854823967 * Q(0)+0.4955056906895569 * Q(1)-0.0763979370214948 * Q(2)+0.0000000000000000 * Q(3)
W4 = +0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+0.3173172288805612 * Q(0)+1.0226557255923103 * Q(1)-0.4278312936233303 * Q(2)+0.0878583391504589 * Q(3)

  W(4)   = (omega1*W1 + omega2*W2 + omega3*W3 +omega4*W4)
!-------------------------------------------------------------------------------!
END SUBROUTINE WENO7_SecondSweep
!===============================================================================!
!
!
!
!
!===============================================================================!
FUNCTION First_Derivative_Central_Order2(f,dx)
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: f(-1:1)
REAL,INTENT(IN) :: dx
REAL            :: First_Derivative_Central_Order2
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!

First_Derivative_Central_Order2 = (f(1)-f(-1))/(2.0*dx)

!-------------------------------------------------------------------------------!
END FUNCTION First_Derivative_Central_Order2
!===============================================================================!
!
!
!
!
!===============================================================================!
FUNCTION Second_Derivative_Central_Order2(f,dx)
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: f(-1:1)
REAL,INTENT(IN) :: dx
REAL            :: Second_Derivative_Central_Order2
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!

Second_Derivative_Central_Order2 = (f(1)-2.0*f(0)+f(-1))/(dx**2)

!-------------------------------------------------------------------------------!
END FUNCTION Second_Derivative_Central_Order2
!===============================================================================!
!
!
!
!
!===============================================================================!
FUNCTION Second_Mixed_Derivative_Central_Order2(f,dx,dy)
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: f(-1:1,-1:1)
REAL,INTENT(IN) :: dx, dy
REAL            :: Second_Mixed_Derivative_Central_Order2
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!

Second_Mixed_Derivative_Central_Order2 = (f(1,1) - f(1,-1) - f(-1,1) + f(-1,-1))/(4.0*dx*dy)

!-------------------------------------------------------------------------------!
END FUNCTION Second_Mixed_Derivative_Central_Order2
!===============================================================================!
!
!
!
!
!===============================================================================!
FUNCTION Laplacian_Central_Order2(f,dx,dy)
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: f(-1:1,-1:1)
REAL,INTENT(IN) :: dx, dy
REAL            :: Laplacian_Central_Order2
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!

Laplacian_Central_Order2 = Second_Derivative_Central_Order2(f(-1:1,0),dx) + Second_Derivative_Central_Order2(f(0,-1:1),dy)

!-------------------------------------------------------------------------------!
END FUNCTION Laplacian_Central_Order2
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ReconstructionX_INPUT_PRIMITIVE(W,vbXs,vbXe,vbYs,vbYe,lbXs,lbXe,lbYs,lbYe,ReconstructionInput,Troubled_Cell_W)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructedVariable
USE MOD_FiniteVolume2D_vars,ONLY: LLX
USE MOD_FiniteVolume2D_vars,ONLY: RRX
USE MOD_FiniteVolume2D_vars,ONLY: LLY
USE MOD_FiniteVolume2D_vars,ONLY: RRY
USE MOD_FiniteVolume2D_vars,ONLY: MCtoP
USE MOD_FiniteVolume2D_vars,ONLY: MPtoC
USE MOD_FiniteVolume2D_vars,ONLY: Tau_Standard_Limiting
USE MOD_FiniteVolume2D_vars,ONLY: Theta_Standard_Limiting
USE MOD_FiniteVolume2D_vars,ONLY: Tau_Overcompressive
USE MOD_FiniteVolume2D_vars,ONLY: Theta_Overcompressive
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,vbXs:vbXe,vbYs:vbYe),    INTENT(IN)          :: W
INTEGER,                                        INTENT(IN)          :: vbXs,vbXe,vbYs,vbYe !*Bounds for the vector
INTEGER,                                        INTENT(IN)          :: lbXs,lbXe,lbYs,lbYe !*Bounds for the loop   
INTEGER, DIMENSION(vbXs:vbXe,vbYs:vbYe),        INTENT(IN),OPTIONAL :: Troubled_Cell_W
INTEGER,                                        INTENT(IN),OPTIONAL :: ReconstructionInput
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj, kk, iGP, indi, indj
INTEGER            :: ReconstructionUsed
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: Wsupport(1:nVar,-nGhosts:nGhosts)
REAL               :: WsupportBox(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts)
REAL               :: LL(1:nVar,1:nVar)
REAL               :: RR(1:nVar,1:nVar)
REAL               :: MM(1:nVar,1:nVar)
INTEGER            :: i_interface
#ifdef MULTIFLUID
REAL               ::   roM,    roP,   roHAT
REAL               ::    uM,     uP,    uHAT
REAL               ::  phiM,   phiP,  phiHAT
REAL               :: pinfM,  pinfP, pinfHAT
REAL               ::    pM,     pP,    pHAT
REAL               ::  gmmM,   gmmP,  gmmHAT
REAL               ::    cM,     cP,    cHAT
REAL               :: Wdummy(1:nVar,1:nGPs)
#endif
!-------------------------------------------------------------------------------!

IF (PRESENT(ReconstructionInput)) THEN
  ReconstructionUsed=ReconstructionInput
ELSE
  ReconstructionUsed=Reconstruction
END IF

SELECT CASE (ReconstructedVariable)
  CASE(0,1,2)

    DO jj=lbYs,lbYe
      DO ii=lbXs,lbXe

        SELECT CASE (ReconstructedVariable)
          CASE(0) !*Conservative
            DO kk=-nGhosts,nGhosts
              Wsupport(1:nVar,kk)=MATMUL(MPtoC(:,:,ii,jj),W(1:nVar,ii+kk,jj)) !*Prim->Cons
            END DO
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=MATMUL(MPtoC(:,:,ii,jj),W(1:nVar,ii+indi,jj+indj))  !*Prim->Cons
              END DO
            END DO
          CASE(1) !*Characteristic
            DO kk=-nGhosts,nGhosts
              Wsupport(1:nVar,kk)=MATMUL( LLX(:,:,ii,jj),W(1:nVar,ii+kk,jj) ) !*Prim->Char
            END DO
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=W(1:nVar,ii+indi,jj+indj)  !*Prim
              END DO
            END DO
          CASE(2) !*Primitive
            Wsupport(1:nVar,-nGhosts:nGhosts)=W(1:nVar,-nGhosts+ii:ii+nGhosts,jj)
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=W(1:nVar,ii+indi,jj+indj)
              END DO
            END DO
          CASE DEFAULT
            ErrorMessage = "Error in ReconstructionX. ReconstructedVariable not available"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT


        SELECT CASE (ABS(ReconstructionUsed))
          CASE(1,10)
            WM(1:nVar,nGPs,ii,jj) = Wsupport(1:nVar,0)
            WP(1:nVar,nGPs,ii,jj) = Wsupport(1:nVar,0)
          CASE(2,20) !*MUSCL
            CALL MUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(21) !*k-MUSCL
            CALL kMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(22) !*MUSCL_CO
            CALL COMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(23) !*MUSCL_VL
            CALL VLMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(24) !*MUSCL_M
            CALL M2MUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(25) !*MUSCL_VA
            CALL VAMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(26) !*MUSCL_SBM
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(27) !*MUSCL_THETA
            CALL THETAMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(28) !*MUSCL_SBM standard
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1),tau_input=Tau_Standard_Limiting,theta_input=Theta_Standard_Limiting)
          CASE(29) !*MUSCL_SBM overcompressive
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1),tau_input=Tau_Overcompressive,theta_input=Theta_Overcompressive)
          CASE(3,4,5,7)
            !*NB: Here I'm entering with primitive variables
            CALL WENO_XDIR(&
                        WsupportBox(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts),& !*Entering with primitive
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(1),&
                        ReconstructionUsed,&
                        RRX(1:nVar,1:nVar,ii,jj),&
                        LLX(1:nVar,1:nVar,ii,jj),&
                        RRY(1:nVar,1:nVar,ii,jj),&
                        LLY(1:nVar,1:nVar,ii,jj))
          CASE DEFAULT
            ErrorMessage = "ReconstructionX_INPUT not implemented"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT


      !*We transform back
        SELECT CASE (ReconstructedVariable)
          CASE(0) !*Conservative

            DO iGP=1,nGPs
              WM(1:nVar,iGP,ii,jj)=MATMUL(MCtoP(:,:,ii,jj),WM(1:nVar,iGP,ii,jj))
              WP(1:nVar,iGP,ii,jj)=MATMUL(MCtoP(:,:,ii,jj),WP(1:nVar,iGP,ii,jj))
            END DO

          CASE(1) !*Characteristic

            SELECT CASE (ABS(ReconstructionUsed))
              CASE(1,10,2,20,21,22,23,24,25,26,27,28,29)

                DO iGP=1,nGPs
                  WM(1:nVar,iGP,ii,jj)=MATMUL(RRX(:,:,ii,jj),WM(1:nVar,iGP,ii,jj))
                  WP(1:nVar,iGP,ii,jj)=MATMUL(RRX(:,:,ii,jj),WP(1:nVar,iGP,ii,jj))
                END DO

              CASE(3,4,5,7)

                !*DO iGP=1,nGPs
                !*  WM(1:nVar,iGP,ii,jj)=WM(1:nVar,iGP,ii,jj)
                !*  WP(1:nVar,iGP,ii,jj)=WP(1:nVar,iGP,ii,jj)
                !*END DO

              CASE DEFAULT

                PRINT*, "YOU SHOULD NOT BE HERE IN ReconstructionX_INPUT_PRIMITIVE"
                STOP

            END SELECT


          CASE(2) !*Primitive

            !*Do nothing, primitive variables have been reconstructed

          CASE DEFAULT
            ErrorMessage = "Error in ReconstructionX_INPUT_PRIMITIVE. ReconstructedVariable not available"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT




#ifdef PWCINTROUBLEDCELLS
#ifdef MULTIFLUID
        IF(PRESENT(Troubled_Cell_W)) THEN
          IF (Troubled_Cell_W(ii,jj) .NE. 0) THEN !*REPLACE BY PWC
            WM(1:nVar,nGPs,ii,jj) = W(1:nVar,ii,jj)
            WP(1:nVar,nGPs,ii,jj) = W(1:nVar,ii,jj)
          END IF
        END IF
#else
        PRINT*, "Impossible to use PWCINTROUBLEDCELLS"
        PRINT*, "Disable this flag"
        PRINT*, "You should not be here"
        STOP
#endif
#endif

#ifdef MULTIFLUID
        WM(nVar+1,1,ii,jj)=W(nVar,ii,jj)
        WP(nVar+1,1,ii,jj)=W(nVar,ii,jj)
#endif


  END DO
END DO




  CASE(-1)

#ifndef MULTIFLUID
    PRINT*, "Reconstruction of Characteristic variables not set up for MULTIFLUID not active"
    STOP
#else

#if defined(MOMENTUMINPRIMITIVEVARIABLES) && !defined(RECONSTRUCTIONVELOCITY)
    PRINT*, "Reconstruction of Characteristic variables not set up for MOMENTUMINPRIMITIVEVARIABLES without RECONSTRUCTIONVELOCITY"
    STOP
#endif

    DO jj=lbYs,lbYe
      DO ii=lbXs,lbXe

        !*LEFT INTERFACE
        i_interface=ii-1

        roM=W(1,i_interface,jj)
        roP=W(1,i_interface+1,jj)  
        !*NB: I always enter with the velocity here
        uM =W(2,i_interface,jj)
        uP =W(2,i_interface+1,jj)  
        !*NB: I always enter with rhophi here
        phiM =W(nVar,i_interface,jj)/W(1,i_interface,jj)
        phiP =W(nVar,i_interface+1,jj)/W(1,i_interface+1,jj)   
        GmmM=Get_Gamma(W(nVar,i_interface,jj))
        pinfM=Get_Pressure_Infinity(W(nVar,i_interface,jj))
        GmmP=Get_Gamma(W(nVar,i_interface+1,jj))
        pinfP=Get_Pressure_Infinity(W(nVar,i_interface+1,jj))
        pM=W(4,i_interface,jj)
        pP=W(4,i_interface+1,jj)  

        !*Roe averages 
        roHAT=SQRT(roM*roP) 
        uHAT =(SQRT(roM)*uM+SQRT(roP)*uP)/(SQRT(roM)+SQRT(roP))                   
        phiHAT=(SQRT(roM)*phiM+SQRT(roP)*phiP)/(SQRT(roM)+SQRT(roP))                    
        pinfHAT=(SQRT(roM)*pinfM+SQRT(roP)*pinfP)/(SQRT(roM)+SQRT(roP))                     
        pHat=(SQRT(roM)*pM+SQRT(roP)*pP)/(SQRT(roM)+SQRT(roP))                     
        gmmHAT=(SQRT(roM)*gmmM+SQRT(roP)*gmmP)/(SQRT(roM)+SQRT(roP))                   
        cHat=SQRT(ABS(gmmHAT*(pHat+pinfHAT)/roHAT)) 

        LL=0.0
        RR=0.0

        RR(1,1)=1.0/cHAT**2
        RR(2,1)=-1.0/(roHat*cHAT)
        RR(4,1)=1.0

        RR(nVar,2)=1.0

        RR(3,3)=1.0

        RR(1,4)=1.0

        RR(1,nVar)=1.0/cHAT**2
        RR(2,nVar)=1.0/(roHat*cHAT)
        RR(4,nVar)=1.0


        LL(4,1)=1.0

        LL(1,2)=-roHAT*cHAT/2.0
        LL(nVar,2)=roHAT*cHAT/2.0

        LL(3,3)=1.0

        LL(1,4)=1.0/2.0
        LL(4,4)=-1.0/cHAT**2
        LL(nVar,4)=1.0/2.0

        LL(2,nVar)=1.0

        DO kk=-nGhosts,nGhosts
          !*NB: I need to have phi in the last component
          Wsupport(1:nVar,kk)=W(1:nVar,ii+kk,jj)
          Wsupport(nVar,kk)=Wsupport(nVar,kk)/Wsupport(1,kk)
          Wsupport(1:nVar,kk)=MATMUL(LL(:,:),Wsupport(1:nVar,kk))
        END DO

        ! MM=MATMUL(LL,RR)
        ! DO indi=1,nVar
        !   PRINT*, indi, MM(indi,:)
        ! END DO
        ! PRINT*


        SELECT CASE (ABS(ReconstructionUsed))
          CASE(1,10)
            WM(1:nVar,nGPs,ii,jj) = Wsupport(1:nVar,0)
          CASE(2,20) !*MUSCL
            CALL MUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      Wdummy(1:nVar,1:nGPs),&
                      MESH_DX(1))
          CASE(21) !*k-MUSCL
            CALL kMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      Wdummy(1:nVar,1:nGPs),&
                      MESH_DX(1))
          CASE(22) !*MUSCL_CO
            CALL COMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      Wdummy(1:nVar,1:nGPs),&
                      MESH_DX(1))
          CASE(23) !*MUSCL_VL
            CALL VLMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      Wdummy(1:nVar,1:nGPs),&
                      MESH_DX(1))
          CASE(24) !*MUSCL_M
            CALL M2MUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      Wdummy(1:nVar,1:nGPs),&
                      MESH_DX(1))
          CASE(25) !*MUSCL_VA
            CALL VAMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      Wdummy(1:nVar,1:nGPs),&
                      MESH_DX(1))
          CASE(26) !*MUSCL_SBM
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      Wdummy(1:nVar,1:nGPs),&
                      MESH_DX(1))
          CASE(27) !*MUSCL_THETA
            CALL THETAMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      Wdummy(1:nVar,1:nGPs),&
                      MESH_DX(1))
          CASE(28) !*MUSCL_SBM standard
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      Wdummy(1:nVar,1:nGPs),&
                      MESH_DX(1),tau_input=Tau_Standard_Limiting,theta_input=Theta_Standard_Limiting)
          CASE(29) !*MUSCL_SBM overcompressive
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      Wdummy(1:nVar,1:nGPs),&
                      MESH_DX(1),tau_input=Tau_Overcompressive,theta_input=Theta_Overcompressive)
          CASE(3,4,5,7)
            !*NB: No shock indicator
            PRINT*, "WENO reconstruction not allowed when reconstructing Chracteristic variables in Active Flux"
            STOP
          CASE DEFAULT
            ErrorMessage = "ReconstructionX_INPUT not implemented"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT

        DO iGP=1,nGPs
          WM(1:nVar,iGP,ii,jj)=MATMUL(RR(:,:),WM(1:nVar,iGP,ii,jj))
          WM(nVar,iGP,ii,jj)=WM(nVar,iGP,ii,jj)*WM(1,iGP,ii,jj)
        END DO

        !*RIGHT INTERFACE
        i_interface=ii

        roM=W(1,i_interface,jj)
        roP=W(1,i_interface+1,jj)  

        !*NB: I always enter with the velocity here
        uM =W(2,i_interface,jj)
        uP =W(2,i_interface+1,jj)  

        phiM =W(nVar,i_interface,jj)/W(1,i_interface,jj)
        phiP =W(nVar,i_interface+1,jj)/W(1,i_interface+1,jj)   
        GmmM=Get_Gamma(W(nVar,i_interface,jj))
        pinfM=Get_Pressure_Infinity(W(nVar,i_interface,jj))
        GmmP=Get_Gamma(W(nVar,i_interface+1,jj))
        pinfP=Get_Pressure_Infinity(W(nVar,i_interface+1,jj))
        pM=W(4,i_interface,jj)
        pP=W(4,i_interface+1,jj)  

        !*Roe averages 
        roHAT=SQRT(roM*roP) 
        uHAT =(SQRT(roM)*uM+SQRT(roP)*uP)/(SQRT(roM)+SQRT(roP))                   
        phiHAT=(SQRT(roM)*phiM+SQRT(roP)*phiP)/(SQRT(roM)+SQRT(roP))                    
        pinfHAT=(SQRT(roM)*pinfM+SQRT(roP)*pinfP)/(SQRT(roM)+SQRT(roP))                     
        pHat=(SQRT(roM)*pM+SQRT(roP)*pP)/(SQRT(roM)+SQRT(roP))                     
        gmmHAT=(SQRT(roM)*gmmM+SQRT(roP)*gmmP)/(SQRT(roM)+SQRT(roP))                   
        cHat=SQRT(ABS(gmmHAT*(pHat+pinfHAT)/roHAT)) 

        LL=0.0
        RR=0.0

        RR(1,1)=1.0/cHAT**2
        RR(2,1)=-1.0/(roHat*cHAT)
        RR(4,1)=1.0

        RR(nVar,2)=1.0

        RR(3,3)=1.0

        RR(1,4)=1.0

        RR(1,nVar)=1.0/cHAT**2
        RR(2,nVar)=1.0/(roHat*cHAT)
        RR(4,nVar)=1.0


        LL(4,1)=1.0

        LL(1,2)=-roHAT*cHAT/2.0
        LL(nVar,2)=roHAT*cHAT/2.0

        LL(3,3)=1.0

        LL(1,4)=1.0/2.0
        LL(4,4)=-1.0/cHAT**2
        LL(nVar,4)=1.0/2.0

        LL(2,nVar)=1.0

        DO kk=-nGhosts,nGhosts
          !*NB: I need to have phi in the last component
          Wsupport(1:nVar,kk)=W(1:nVar,ii+kk,jj)
          Wsupport(nVar,kk)=Wsupport(nVar,kk)/Wsupport(1,kk)
          Wsupport(1:nVar,kk)=MATMUL(LL(:,:),Wsupport(1:nVar,kk))
        END DO


        SELECT CASE (ABS(ReconstructionUsed))
          CASE(1,10)
            WP(1:nVar,nGPs,ii,jj) = Wsupport(1:nVar,0)
          CASE(2,20) !*MUSCL
            CALL MUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      Wdummy(1:nVar,1:nGPs),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(21) !*k-MUSCL
            CALL kMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      Wdummy(1:nVar,1:nGPs),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(22) !*MUSCL_CO
            CALL COMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      Wdummy(1:nVar,1:nGPs),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(23) !*MUSCL_VL
            CALL VLMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      Wdummy(1:nVar,1:nGPs),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(24) !*MUSCL_M
            CALL M2MUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      Wdummy(1:nVar,1:nGPs),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(25) !*MUSCL_VA
            CALL VAMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      Wdummy(1:nVar,1:nGPs),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(26) !*MUSCL_SBM
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      Wdummy(1:nVar,1:nGPs),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(27) !*MUSCL_THETA
            CALL THETAMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      Wdummy(1:nVar,1:nGPs),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(28) !*MUSCL_SBM standard
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      Wdummy(1:nVar,1:nGPs),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1),tau_input=Tau_Standard_Limiting,theta_input=Theta_Standard_Limiting)
          CASE(29) !*MUSCL_SBM overcompressive
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      Wdummy(1:nVar,1:nGPs),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1),tau_input=Tau_Overcompressive,theta_input=Theta_Overcompressive)
          CASE(3,4,5,7)
            !*NB: No shock indicator
            PRINT*, "WENO reconstruction not allowed when reconstructing Chracteristic variables in Active Flux"
            STOP
          CASE DEFAULT
            ErrorMessage = "ReconstructionX_INPUT not implemented"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT


        DO iGP=1,nGPs
          WP(1:nVar,iGP,ii,jj)=MATMUL(RR(:,:),WP(1:nVar,iGP,ii,jj))
          WP(nVar,iGP,ii,jj)=WP(nVar,iGP,ii,jj)*WP(1,iGP,ii,jj)
        END DO



#ifdef PWCINTROUBLEDCELLS
#ifdef MULTIFLUID
        IF(PRESENT(Troubled_Cell_W)) THEN
          IF (Troubled_Cell_W(ii,jj) .NE. 0) THEN !*REPLACE BY PWC
            WM(1:nVar,nGPs,ii,jj) = W(1:nVar,ii,jj)
            WP(1:nVar,nGPs,ii,jj) = W(1:nVar,ii,jj)
          END IF
        END IF
#else
        PRINT*, "Impossible to use PWCINTROUBLEDCELLS"
        PRINT*, "Disable this flag"
        PRINT*, "You should not be here"
        STOP
#endif
#endif


#ifdef MULTIFLUID
        WM(nVar+1,1,ii,jj)=W(nVar,ii,jj)
        WP(nVar+1,1,ii,jj)=W(nVar,ii,jj)
#endif


  END DO
END DO


#endif

  CASE DEFAULT
    ErrorMessage = "Error in ReconstructionX_INPUT. ReconstructedVariable not available"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT




!-------------------------------------------------------------------------------!
END SUBROUTINE ReconstructionX_INPUT_PRIMITIVE
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ReconstructionY_INPUT_PRIMITIVE(W,vbXs,vbXe,vbYs,vbYe,lbXs,lbXe,lbYs,lbYe,ReconstructionInput,Troubled_Cell_W)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructedVariable
USE MOD_FiniteVolume2D_vars,ONLY: LLX
USE MOD_FiniteVolume2D_vars,ONLY: RRX
USE MOD_FiniteVolume2D_vars,ONLY: LLY
USE MOD_FiniteVolume2D_vars,ONLY: RRY
USE MOD_FiniteVolume2D_vars,ONLY: MCtoP
USE MOD_FiniteVolume2D_vars,ONLY: MPtoC
USE MOD_FiniteVolume2D_vars,ONLY: Tau_Standard_Limiting
USE MOD_FiniteVolume2D_vars,ONLY: Theta_Standard_Limiting
USE MOD_FiniteVolume2D_vars,ONLY: Tau_Overcompressive
USE MOD_FiniteVolume2D_vars,ONLY: Theta_Overcompressive
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,vbXs:vbXe,vbYs:vbYe),    INTENT(IN)          :: W
INTEGER,                                        INTENT(IN)          :: vbXs,vbXe,vbYs,vbYe !*Bounds for the vector
INTEGER,                                        INTENT(IN)          :: lbXs,lbXe,lbYs,lbYe !*Bounds for the loop    
INTEGER, DIMENSION(vbXs:vbXe,vbYs:vbYe),        INTENT(IN),OPTIONAL :: Troubled_Cell_W                        
INTEGER,                                        INTENT(IN),OPTIONAL :: ReconstructionInput
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj, kk, iGP
INTEGER            :: ReconstructionUsed
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: Wsupport(1:nVar,-nGhosts:nGhosts)
REAL               :: WsupportBox(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts)
REAL               :: LL(1:nVar,1:nVar)
REAL               :: RR(1:nVar,1:nVar)
REAL               :: MM(1:nVar,1:nVar)
INTEGER            :: j_interface
REAL               :: a_support, b_support
#ifdef MULTIFLUID
REAL               ::   roM,    roP,   roHAT
REAL               ::    uM,     uP,    uHAT
REAL               ::  phiM,   phiP,  phiHAT
REAL               :: pinfM,  pinfP, pinfHAT
REAL               ::    pM,     pP,    pHAT
REAL               ::  gmmM,   gmmP,  gmmHAT
REAL               ::    cM,     cP,    cHAT
REAL               :: Wdummy(1:nVar,1:nGPs)
#endif
INTEGER            :: indi, indj
!-------------------------------------------------------------------------------!

IF (PRESENT(ReconstructionInput)) THEN
  ReconstructionUsed=ReconstructionInput
ELSE
  ReconstructionUsed=Reconstruction
END IF


SELECT CASE (ReconstructedVariable)
  CASE(0,1,2)

    DO jj=lbYs,lbYe
      DO ii=lbXs,lbXe

        SELECT CASE (ReconstructedVariable)
          CASE(0) !*Conservative
            DO kk=-nGhosts,nGhosts
              Wsupport(1:nVar,kk)=MATMUL(MPtoC(:,:,ii,jj),W(1:nVar,ii,jj+kk)) !*Prim->Cons
            END DO
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=MATMUL(MPtoC(:,:,ii,jj),W(1:nVar,ii+indi,jj+indj))  !*Prim->Cons
              END DO
            END DO
          CASE(1) !*Characteristic
            DO kk=-nGhosts,nGhosts
              Wsupport(1:nVar,kk)=MATMUL( LLY(:,:,ii,jj), W(1:nVar,ii,jj+kk) ) !*Prim->Char
            END DO
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=W(1:nVar,ii+indi,jj+indj)  !*Prim
              END DO
            END DO
          CASE(2) !*Primitive
            Wsupport(1:nVar,-nGhosts:nGhosts)=W(1:nVar,ii,-nGhosts+jj:jj+nGhosts)
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=W(1:nVar,ii+indi,jj+indj)
              END DO
            END DO
          CASE DEFAULT
            ErrorMessage = "Error in ReconstructionY. ReconstructedVariable not available"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT



        SELECT CASE (ABS(ReconstructionUsed))
          CASE(1,10)
            WM(1:nVar,nGPs,ii,jj) = Wsupport(1:nVar,0)
            WP(1:nVar,nGPs,ii,jj) = Wsupport(1:nVar,0)
          CASE(2,20) !*MUSCL
            CALL MUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(21) !*k-MUSCL
            CALL kMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(22) !*MUSCL_CO
            CALL COMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(23) !*MUSCL_VL
            CALL VLMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(24) !*MUSCL_M
            CALL M2MUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(25) !*MUSCL_VA
            CALL VAMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(26) !*MUSCL_SBM
            CALL SBMMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(27) !*MUSCL_THETA
            CALL THETAMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(28) !*MUSCL_SBM standard
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(2),tau_input=Tau_Standard_Limiting,theta_input=Theta_Standard_Limiting)
          CASE(29) !*MUSCL_SBM overcompressive
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(2),tau_input=Tau_Overcompressive,theta_input=Theta_Overcompressive)
          CASE(3,4,5,7)
            !*NB: Here I'm entering with primitive variables
            CALL WENO_YDIR(&
                      WsupportBox(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts),& !*Entering with Prim
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(2),&
                      ReconstructionUsed,&
                      RRX(1:nVar,1:nVar,ii,jj),&
                      LLX(1:nVar,1:nVar,ii,jj),&
                      RRY(1:nVar,1:nVar,ii,jj),&
                      LLY(1:nVar,1:nVar,ii,jj))
          CASE DEFAULT
            ErrorMessage = "Reconstruction not implemented in ReconstructionY_INPUT"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT

      !*We transform back
        SELECT CASE (ReconstructedVariable)
          CASE(0) !*Conservative

            DO iGP=1,nGPs
              WM(1:nVar,iGP,ii,jj)=MATMUL(MCtoP(:,:,ii,jj),WM(1:nVar,iGP,ii,jj))
              WP(1:nVar,iGP,ii,jj)=MATMUL(MCtoP(:,:,ii,jj),WP(1:nVar,iGP,ii,jj))
            END DO

          CASE(1) !*Characteristic

            SELECT CASE (ABS(ReconstructionUsed))
              CASE(1,10,2,20,21,22,23,24,25,26,27,28,29)

                DO iGP=1,nGPs
                  WM(1:nVar,iGP,ii,jj)=MATMUL(RRY(:,:,ii,jj),WM(1:nVar,iGP,ii,jj))
                  WP(1:nVar,iGP,ii,jj)=MATMUL(RRY(:,:,ii,jj),WP(1:nVar,iGP,ii,jj))
                END DO

              CASE(3,4,5,7)

                !*DO iGP=1,nGPs
                !*  WM(1:nVar,iGP,ii,jj)=WM(1:nVar,iGP,ii,jj)
                !*  WP(1:nVar,iGP,ii,jj)=WP(1:nVar,iGP,ii,jj)
                !*END DO

              CASE DEFAULT

                PRINT*, "YOU SHOULD NOT BE HERE IN ReconstructionY_INPUT_PRIMITIVE"
                STOP

            END SELECT


          CASE(2) !*Primitive

            !*Do nothing, primitive variables have been reconstructed

          CASE DEFAULT
            ErrorMessage = "Error in ReconstructionY_INPUT_PRIMITIVE. ReconstructedVariable not available"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT





#ifdef PWCINTROUBLEDCELLS
#ifdef MULTIFLUID
        IF(PRESENT(Troubled_Cell_W)) THEN
          IF (Troubled_Cell_W(ii,jj) .NE. 0) THEN !*REPLACE BY PWC
            WM(1:nVar,nGPs,ii,jj) = W(1:nVar,ii,jj)
            WP(1:nVar,nGPs,ii,jj) = W(1:nVar,ii,jj)
          END IF
        END IF
#else
        PRINT*, "Impossible to use PWCINTROUBLEDCELLS"
        PRINT*, "Disable this flag"
        PRINT*, "You should not be here"
        STOP
#endif
#endif

#ifdef MULTIFLUID
        WM(nVar+1,1,ii,jj)=W(nVar,ii,jj)
        WP(nVar+1,1,ii,jj)=W(nVar,ii,jj)
#endif


      END DO
    END DO

  CASE(-1)

#ifndef MULTIFLUID
    PRINT*, "Reconstruction of Characteristic variables not set up for MULTIFLUID not active"
    STOP
#else

#if defined(MOMENTUMINPRIMITIVEVARIABLES) && !defined(RECONSTRUCTIONVELOCITY)
    PRINT*, "Reconstruction of Characteristic variables not set up for MOMENTUMINPRIMITIVEVARIABLES without RECONSTRUCTIONVELOCITY"
    STOP
#endif

    DO jj=lbYs,lbYe
      DO ii=lbXs,lbXe


        !*BOTTOM INTERFACE
        j_interface=jj-1

        roM=W(1,ii,j_interface)
        roP=W(1,ii,j_interface+1)  
        !*NB: I always enter with the velocity here
        uM =W(3,ii,j_interface)   !*NB: This time I take the v component
        uP =W(3,ii,j_interface+1) !*NB: This time I take the v component 
        !*NB: I always enter with rhophi here
        phiM =W(nVar,ii,j_interface)/W(1,ii,j_interface)
        phiP =W(nVar,ii,j_interface+1)/W(1,ii,j_interface+1)   
        GmmM=Get_Gamma(W(nVar,ii,j_interface))
        pinfM=Get_Pressure_Infinity(W(nVar,ii,j_interface))
        GmmP=Get_Gamma(W(nVar,ii,j_interface+1))
        pinfP=Get_Pressure_Infinity(W(nVar,ii,j_interface+1))
        pM=W(4,ii,j_interface)
        pP=W(4,ii,j_interface+1)  

        !*Roe averages 
        roHAT=SQRT(roM*roP) 
        uHAT =(SQRT(roM)*uM+SQRT(roP)*uP)/(SQRT(roM)+SQRT(roP))                   
        phiHAT=(SQRT(roM)*phiM+SQRT(roP)*phiP)/(SQRT(roM)+SQRT(roP))                    
        pinfHAT=(SQRT(roM)*pinfM+SQRT(roP)*pinfP)/(SQRT(roM)+SQRT(roP))                     
        pHat=(SQRT(roM)*pM+SQRT(roP)*pP)/(SQRT(roM)+SQRT(roP))                     
        gmmHAT=(SQRT(roM)*gmmM+SQRT(roP)*gmmP)/(SQRT(roM)+SQRT(roP))                   
        cHat=SQRT(ABS(gmmHAT*(pHat+pinfHAT)/roHAT)) 

        LL=0.0
        RR=0.0

        RR(1,1)=1.0/cHAT**2
        RR(2,1)=-1.0/(roHat*cHAT)
        RR(4,1)=1.0

        RR(nVar,2)=1.0

        RR(3,3)=1.0

        RR(1,4)=1.0

        RR(1,nVar)=1.0/cHAT**2
        RR(2,nVar)=1.0/(roHat*cHAT)
        RR(4,nVar)=1.0


        LL(4,1)=1.0

        LL(1,2)=-roHAT*cHAT/2.0
        LL(nVar,2)=roHAT*cHAT/2.0

        LL(3,3)=1.0

        LL(1,4)=1.0/2.0
        LL(4,4)=-1.0/cHAT**2
        LL(nVar,4)=1.0/2.0

        LL(2,nVar)=1.0

        DO kk=-nGhosts,nGhosts
          Wsupport(1:nVar,kk)=W(1:nVar,ii,jj+kk)
          !*NB: Switching velocities
          Wsupport(2,kk)=W(3,ii,jj+kk) 
          Wsupport(3,kk)=W(2,ii,jj+kk) 
          !*NB: I need to have phi in the last component
          Wsupport(nVar,kk)=Wsupport(nVar,kk)/Wsupport(1,kk)
          Wsupport(1:nVar,kk)=MATMUL(LL(:,:),Wsupport(1:nVar,kk))
        END DO



        SELECT CASE (ABS(ReconstructionUsed))
          CASE(1,10)
            WM(1:nVar,nGPs,ii,jj) = Wsupport(1:nVar,0)
          CASE(2,20) !*MUSCL
            CALL MUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        Wdummy(1:nVar,1:nGPs),&
                        MESH_DX(2))
          CASE(21) !*k-MUSCL
            CALL kMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        Wdummy(1:nVar,1:nGPs),&
                        MESH_DX(2))
          CASE(22) !*MUSCL_CO
            CALL COMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        Wdummy(1:nVar,1:nGPs),&
                        MESH_DX(2))
          CASE(23) !*MUSCL_VL
            CALL VLMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        Wdummy(1:nVar,1:nGPs),&
                        MESH_DX(2))
          CASE(24) !*MUSCL_M
            CALL M2MUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        Wdummy(1:nVar,1:nGPs),&
                        MESH_DX(2))
          CASE(25) !*MUSCL_VA
            CALL VAMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        Wdummy(1:nVar,1:nGPs),&
                        MESH_DX(2))
          CASE(26) !*MUSCL_SBM
            CALL SBMMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        Wdummy(1:nVar,1:nGPs),&
                        MESH_DX(2))
          CASE(27) !*MUSCL_THETA
            CALL THETAMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        Wdummy(1:nVar,1:nGPs),&
                        MESH_DX(2))
          CASE(28) !*MUSCL_SBM standard
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      Wdummy(1:nVar,1:nGPs),&
                      MESH_DX(2),tau_input=Tau_Standard_Limiting,theta_input=Theta_Standard_Limiting)
          CASE(29) !*MUSCL_SBM overcompressive
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      Wdummy(1:nVar,1:nGPs),&
                      MESH_DX(2),tau_input=Tau_Overcompressive,theta_input=Theta_Overcompressive)
          CASE(3,4,5,7)
            !*NB: No shock indicator
            PRINT*, "WENO reconstruction not allowed when reconstructing Chracteristic variables in Active Flux"
            STOP
          CASE DEFAULT
            ErrorMessage = "Reconstruction not implemented in ReconstructionY_INPUT"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT

        DO iGP=1,nGPs
          WM(1:nVar,iGP,ii,jj)=MATMUL(RR(:,:),WM(1:nVar,iGP,ii,jj))
          WM(nVar,iGP,ii,jj)=WM(nVar,iGP,ii,jj)*WM(1,iGP,ii,jj)
          !*NB: Switching back the velocities
          a_support=WM(2,iGP,ii,jj) 
          b_support=WM(3,iGP,ii,jj) 
          WM(2,iGP,ii,jj)=b_support 
          WM(3,iGP,ii,jj)=a_support 
        END DO

        !*UPPER INTERFACE
        j_interface=jj

        roM=W(1,ii,j_interface)
        roP=W(1,ii,j_interface+1)  
        !*NB: I always enter with the velocity here
        uM =W(3,ii,j_interface)   !*NB: This time I take the v component
        uP =W(3,ii,j_interface+1) !*NB: This time I take the v component 
        !*NB: I always enter with rhophi here
        phiM =W(nVar,ii,j_interface)/W(1,ii,j_interface)
        phiP =W(nVar,ii,j_interface+1)/W(1,ii,j_interface+1)   
        GmmM=Get_Gamma(W(nVar,ii,j_interface))
        pinfM=Get_Pressure_Infinity(W(nVar,ii,j_interface))
        GmmP=Get_Gamma(W(nVar,ii,j_interface+1))
        pinfP=Get_Pressure_Infinity(W(nVar,ii,j_interface+1))
        pM=W(4,ii,j_interface)
        pP=W(4,ii,j_interface+1)  

        !*Roe averages 
        roHAT=SQRT(roM*roP) 
        uHAT =(SQRT(roM)*uM+SQRT(roP)*uP)/(SQRT(roM)+SQRT(roP))                   
        phiHAT=(SQRT(roM)*phiM+SQRT(roP)*phiP)/(SQRT(roM)+SQRT(roP))                    
        pinfHAT=(SQRT(roM)*pinfM+SQRT(roP)*pinfP)/(SQRT(roM)+SQRT(roP))                     
        pHat=(SQRT(roM)*pM+SQRT(roP)*pP)/(SQRT(roM)+SQRT(roP))                     
        gmmHAT=(SQRT(roM)*gmmM+SQRT(roP)*gmmP)/(SQRT(roM)+SQRT(roP))                   
        cHat=SQRT(ABS(gmmHAT*(pHat+pinfHAT)/roHAT)) 

        LL=0.0
        RR=0.0

        RR(1,1)=1.0/cHAT**2
        RR(2,1)=-1.0/(roHat*cHAT)
        RR(4,1)=1.0

        RR(nVar,2)=1.0

        RR(3,3)=1.0

        RR(1,4)=1.0

        RR(1,nVar)=1.0/cHAT**2
        RR(2,nVar)=1.0/(roHat*cHAT)
        RR(4,nVar)=1.0


        LL(4,1)=1.0

        LL(1,2)=-roHAT*cHAT/2.0
        LL(nVar,2)=roHAT*cHAT/2.0

        LL(3,3)=1.0

        LL(1,4)=1.0/2.0
        LL(4,4)=-1.0/cHAT**2
        LL(nVar,4)=1.0/2.0

        LL(2,nVar)=1.0

        DO kk=-nGhosts,nGhosts
          Wsupport(1:nVar,kk)=W(1:nVar,ii,jj+kk)
          !*NB: Switching velocities
          Wsupport(2,kk)=W(3,ii,jj+kk) 
          Wsupport(3,kk)=W(2,ii,jj+kk) 
          !*NB: I need to have phi in the last component
          Wsupport(nVar,kk)=Wsupport(nVar,kk)/Wsupport(1,kk)
          Wsupport(1:nVar,kk)=MATMUL(LL(:,:),Wsupport(1:nVar,kk))
        END DO

        SELECT CASE (ABS(ReconstructionUsed))
          CASE(1,10)
            WP(1:nVar,nGPs,ii,jj) = Wsupport(1:nVar,0)
          CASE(2,20) !*MUSCL
            CALL MUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        Wdummy(1:nVar,1:nGPs),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(21) !*k-MUSCL
            CALL kMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        Wdummy(1:nVar,1:nGPs),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(22) !*MUSCL_CO
            CALL COMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        Wdummy(1:nVar,1:nGPs),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(23) !*MUSCL_VL
            CALL VLMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        Wdummy(1:nVar,1:nGPs),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(24) !*MUSCL_M
            CALL M2MUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        Wdummy(1:nVar,1:nGPs),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(25) !*MUSCL_VA
            CALL VAMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        Wdummy(1:nVar,1:nGPs),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(26) !*MUSCL_SBM
            CALL SBMMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        Wdummy(1:nVar,1:nGPs),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(27) !*MUSCL_THETA
            CALL THETAMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        Wdummy(1:nVar,1:nGPs),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(28) !*MUSCL_SBM standard
            CALL SBMMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        Wdummy(1:nVar,1:nGPs),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2),tau_input=Tau_Standard_Limiting,theta_input=Theta_Standard_Limiting)
          CASE(29) !*MUSCL_SBM overcompressive
            CALL SBMMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        Wdummy(1:nVar,1:nGPs),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2),tau_input=Tau_Overcompressive,theta_input=Theta_Overcompressive)
          CASE(3,4,5,7)
            !*NB: No shock indicator
            PRINT*, "WENO reconstruction not allowed when reconstructing Chracteristic variables in Active Flux"
            STOP
          CASE DEFAULT
            ErrorMessage = "Reconstruction not implemented in ReconstructionY_INPUT"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT



        DO iGP=1,nGPs
          WP(1:nVar,iGP,ii,jj)=MATMUL(RR(:,:),WP(1:nVar,iGP,ii,jj))
          WP(nVar,iGP,ii,jj)=WP(nVar,iGP,ii,jj)*WP(1,iGP,ii,jj)
          !*NB: Switching back the velocities
          a_support=WP(2,iGP,ii,jj) 
          b_support=WP(3,iGP,ii,jj) 
          WP(2,iGP,ii,jj)=b_support 
          WP(3,iGP,ii,jj)=a_support 
        END DO



#ifdef PWCINTROUBLEDCELLS
#ifdef MULTIFLUID
        IF(PRESENT(Troubled_Cell_W)) THEN
          IF (Troubled_Cell_W(ii,jj) .NE. 0) THEN !*REPLACE BY PWC
            WM(1:nVar,nGPs,ii,jj) = W(1:nVar,ii,jj)
            WP(1:nVar,nGPs,ii,jj) = W(1:nVar,ii,jj)
          END IF
        END IF
#else
        PRINT*, "Impossible to use PWCINTROUBLEDCELLS"
        PRINT*, "Disable this flag"
        PRINT*, "You should not be here"
        STOP
#endif
#endif

#ifdef MULTIFLUID
        WM(nVar+1,1,ii,jj)=W(nVar,ii,jj)
        WP(nVar+1,1,ii,jj)=W(nVar,ii,jj)
#endif




  END DO
END DO


#endif

  CASE DEFAULT
    ErrorMessage = "Error in ReconstructionY_INPUT. ReconstructedVariable not available"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT




!-------------------------------------------------------------------------------!
END SUBROUTINE ReconstructionY_INPUT_PRIMITIVE
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ReconstructionX_INPUT_CONSERVED(W,vbXs,vbXe,vbYs,vbYe,lbXs,lbXe,lbYs,lbYe,Troubled_Cell_W)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructedVariable
USE MOD_FiniteVolume2D_vars,ONLY: LLX
USE MOD_FiniteVolume2D_vars,ONLY: RRX
USE MOD_FiniteVolume2D_vars,ONLY: LLY
USE MOD_FiniteVolume2D_vars,ONLY: RRY
USE MOD_FiniteVolume2D_vars,ONLY: MCtoP
USE MOD_FiniteVolume2D_vars,ONLY: MPtoC
USE MOD_FiniteVolume2D_vars,ONLY: Tau_Standard_Limiting
USE MOD_FiniteVolume2D_vars,ONLY: Theta_Standard_Limiting
USE MOD_FiniteVolume2D_vars,ONLY: Tau_Overcompressive
USE MOD_FiniteVolume2D_vars,ONLY: Theta_Overcompressive
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,vbXs:vbXe,vbYs:vbYe),    INTENT(IN)          :: W
INTEGER,                                        INTENT(IN)          :: vbXs,vbXe,vbYs,vbYe !*Bounds for the vector
INTEGER,                                        INTENT(IN)          :: lbXs,lbXe,lbYs,lbYe !*Bounds for the loop   
INTEGER, DIMENSION(vbXs:vbXe,vbYs:vbYe),        INTENT(IN),OPTIONAL :: Troubled_Cell_W                        
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj, kk, iGP, indi, indj
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: Wsupport(1:nVar,-nGhosts:nGhosts)
REAL               :: WsupportBox(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts)
REAL               :: LL(1:nVar,1:nVar)
REAL               :: RR(1:nVar,1:nVar)
REAL               :: MM(1:nVar,1:nVar)
INTEGER            :: i_interface
#ifdef MULTIFLUID
REAL               ::   roM,    roP,   roHAT
REAL               ::    uM,     uP,    uHAT
REAL               ::  phiM,   phiP,  phiHAT
REAL               :: pinfM,  pinfP, pinfHAT
REAL               ::    pM,     pP,    pHAT
REAL               ::  gmmM,   gmmP,  gmmHAT
REAL               ::    cM,     cP,    cHAT
REAL               :: Wdummy(1:nVar,1:nGPs)
#endif
!-------------------------------------------------------------------------------!


SELECT CASE (ReconstructedVariable)
  CASE(0,1,2)

    DO jj=lbYs,lbYe
      DO ii=lbXs,lbXe

        SELECT CASE (ReconstructedVariable)
          CASE(0) !*Conserved
            Wsupport(1:nVar,-nGhosts:nGhosts)=W(1:nVar,-nGhosts+ii:ii+nGhosts,jj)
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=W(1:nVar,ii+indi,jj+indj)
              END DO
            END DO
          CASE(1) !*Characteristic
            DO kk=-nGhosts,nGhosts
              Wsupport(1:nVar,kk)=MATMUL( LLX(:,:,ii,jj),W(1:nVar,ii+kk,jj) ) !*Cons->Char
            END DO
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=W(1:nVar,ii+indi,jj+indj)  !*Cons
              END DO
            END DO
          CASE(2) !*Primitive
            DO kk=-nGhosts,nGhosts
              Wsupport(1:nVar,kk)=MATMUL(MCtoP(:,:,ii,jj),W(1:nVar,ii+kk,jj)) !*Cons->Prim
            END DO
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=MATMUL(MCtoP(:,:,ii,jj),W(1:nVar,ii+indi,jj+indj))  !*Cons->Prim
              END DO
            END DO
          CASE DEFAULT
            ErrorMessage = "Error in ReconstructionX. ReconstructedVariable not available"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT


        SELECT CASE (ABS(Reconstruction))
          CASE(1,10)
            WM(1:nVar,nGPs,ii,jj) = Wsupport(1:nVar,0)
            WP(1:nVar,nGPs,ii,jj) = Wsupport(1:nVar,0)
          CASE(2,20) !*MUSCL
            CALL MUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(21) !*k-MUSCL
            CALL kMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(22) !*MUSCL_CO
            CALL COMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(23) !*MUSCL_VL
            CALL VLMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(24) !*MUSCL_M
            CALL M2MUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(25) !*MUSCL_VA
            CALL VAMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(26) !*MUSCL_SBM
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(27) !*MUSCL_THETA
            CALL THETAMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(28) !*MUSCL_SBM standard
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1),tau_input=Tau_Standard_Limiting,theta_input=Theta_Standard_Limiting)
          CASE(29) !*MUSCL_SBM overcompressive
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1),tau_input=Tau_Overcompressive,theta_input=Theta_Overcompressive)
          CASE(3,4,5,7)
            !*NB: Here I'm entering with conserved variables
            CALL WENO_XDIR(&
                        WsupportBox(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts),& !*Entering with Cons
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(1),&
                        Reconstruction,&
                        RRX(1:nVar,1:nVar,ii,jj),&
                        LLX(1:nVar,1:nVar,ii,jj),&
                        RRY(1:nVar,1:nVar,ii,jj),&
                        LLY(1:nVar,1:nVar,ii,jj))
          CASE DEFAULT
            ErrorMessage = "ReconstructionX_INPUT not implemented"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT


      !*We transform back
        SELECT CASE (ReconstructedVariable)
          CASE(0) !*Conservative

            !*Do nothing

          CASE(1) !*Characteristic

            SELECT CASE (ABS(Reconstruction))
              CASE(1,10,2,20,21,22,23,24,25,26,27,28,29)

                DO iGP=1,nGPs
                  WM(1:nVar,iGP,ii,jj)=MATMUL(RRX(:,:,ii,jj),WM(1:nVar,iGP,ii,jj))
                  WP(1:nVar,iGP,ii,jj)=MATMUL(RRX(:,:,ii,jj),WP(1:nVar,iGP,ii,jj))
                END DO

              CASE(3,4,5,7)

                !*Do nothing

              CASE DEFAULT

                PRINT*, "YOU SHOULD NOT BE HERE IN ReconstructionX_INPUT_CONSERVED"
                STOP

            END SELECT


          CASE(2) !*Primitive

                DO iGP=1,nGPs
                  WM(1:nVar,iGP,ii,jj)=MATMUL(MPtoC(:,:,ii,jj),WM(1:nVar,iGP,ii,jj))
                  WP(1:nVar,iGP,ii,jj)=MATMUL(MPtoC(:,:,ii,jj),WP(1:nVar,iGP,ii,jj))
                END DO

          CASE DEFAULT
            ErrorMessage = "Error in ReconstructionX_INPUT_CONSERVED. ReconstructedVariable not available"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT




#ifdef PWCINTROUBLEDCELLS
#ifdef MULTIFLUID
        IF(PRESENT(Troubled_Cell_W)) THEN
          IF (Troubled_Cell_W(ii,jj) .NE. 0) THEN !*REPLACE BY PWC
            WM(1:nVar,nGPs,ii,jj) = W(1:nVar,ii,jj)
            WP(1:nVar,nGPs,ii,jj) = W(1:nVar,ii,jj)
          END IF
        END IF
#else
        PRINT*, "Impossible to use PWCINTROUBLEDCELLS"
        PRINT*, "Disable this flag"
        PRINT*, "You should not be here"
        STOP
#endif
#endif

#ifdef MULTIFLUID
        WM(nVar+1,1,ii,jj)=W(nVar,ii,jj)
        WP(nVar+1,1,ii,jj)=W(nVar,ii,jj)
#endif


  END DO
END DO




  CASE(-1)

    PRINT*, "Reconstruction of Characteristic variables of Alex' and Shaoshuai not set up for subroutine ReconstructionX_INPUT_CONSERVED"
    STOP

  CASE DEFAULT
    ErrorMessage = "Error in ReconstructionX_INPUT. ReconstructedVariable not available"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT




!-------------------------------------------------------------------------------!
END SUBROUTINE ReconstructionX_INPUT_CONSERVED
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ReconstructionY_INPUT_CONSERVED(W,vbXs,vbXe,vbYs,vbYe,lbXs,lbXe,lbYs,lbYe,Troubled_Cell_W)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructedVariable
USE MOD_FiniteVolume2D_vars,ONLY: LLX
USE MOD_FiniteVolume2D_vars,ONLY: RRX
USE MOD_FiniteVolume2D_vars,ONLY: LLY
USE MOD_FiniteVolume2D_vars,ONLY: RRY
USE MOD_FiniteVolume2D_vars,ONLY: MCtoP
USE MOD_FiniteVolume2D_vars,ONLY: MPtoC
USE MOD_FiniteVolume2D_vars,ONLY: Tau_Standard_Limiting
USE MOD_FiniteVolume2D_vars,ONLY: Theta_Standard_Limiting
USE MOD_FiniteVolume2D_vars,ONLY: Tau_Overcompressive
USE MOD_FiniteVolume2D_vars,ONLY: Theta_Overcompressive
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,vbXs:vbXe,vbYs:vbYe),    INTENT(IN)          :: W
INTEGER,                                        INTENT(IN)          :: vbXs,vbXe,vbYs,vbYe !*Bounds for the vector
INTEGER,                                        INTENT(IN)          :: lbXs,lbXe,lbYs,lbYe !*Bounds for the loop    
INTEGER, DIMENSION(vbXs:vbXe,vbYs:vbYe),        INTENT(IN),OPTIONAL :: Troubled_Cell_W                        

!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj, kk, iGP
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: Wsupport(1:nVar,-nGhosts:nGhosts)
REAL               :: WsupportBox(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts)
REAL               :: LL(1:nVar,1:nVar)
REAL               :: RR(1:nVar,1:nVar)
REAL               :: MM(1:nVar,1:nVar)
INTEGER            :: j_interface
REAL               :: a_support, b_support
#ifdef MULTIFLUID
REAL               ::   roM,    roP,   roHAT
REAL               ::    uM,     uP,    uHAT
REAL               ::  phiM,   phiP,  phiHAT
REAL               :: pinfM,  pinfP, pinfHAT
REAL               ::    pM,     pP,    pHAT
REAL               ::  gmmM,   gmmP,  gmmHAT
REAL               ::    cM,     cP,    cHAT
REAL               :: Wdummy(1:nVar,1:nGPs)
#endif
INTEGER            :: indi, indj
!-------------------------------------------------------------------------------!

SELECT CASE (ReconstructedVariable)
  CASE(0,1,2)

    DO jj=lbYs,lbYe
      DO ii=lbXs,lbXe

        SELECT CASE (ReconstructedVariable)
          CASE(0) !*Conserved
            Wsupport(1:nVar,-nGhosts:nGhosts)=W(1:nVar,ii,-nGhosts+jj:jj+nGhosts)
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=W(1:nVar,ii+indi,jj+indj)
              END DO
            END DO

          CASE(1) !*Characteristic
            DO kk=-nGhosts,nGhosts
              Wsupport(1:nVar,kk)=MATMUL( LLY(:,:,ii,jj),W(1:nVar,ii,jj+kk)) !*Cons->Char
            END DO
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=W(1:nVar,ii+indi,jj+indj) !*Cons
              END DO
            END DO

          CASE(2) !*Primitive
            DO kk=-nGhosts,nGhosts
              Wsupport(1:nVar,kk)=MATMUL(MCtoP(:,:,ii,jj),W(1:nVar,ii,jj+kk)) !*Cons->Prim
            END DO
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=MATMUL(MCtoP(:,:,ii,jj),W(1:nVar,ii+indi,jj+indj))  !*Cons->Prim
              END DO
            END DO


          CASE DEFAULT
            ErrorMessage = "Error in ReconstructionY. ReconstructedVariable not available"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT



        SELECT CASE (ABS(Reconstruction))
          CASE(1,10)
            WM(1:nVar,nGPs,ii,jj) = Wsupport(1:nVar,0)
            WP(1:nVar,nGPs,ii,jj) = Wsupport(1:nVar,0)
          CASE(2,20) !*MUSCL
            CALL MUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(21) !*k-MUSCL
            CALL kMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(22) !*MUSCL_CO
            CALL COMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(23) !*MUSCL_VL
            CALL VLMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(24) !*MUSCL_M
            CALL M2MUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(25) !*MUSCL_VA
            CALL VAMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(26) !*MUSCL_SBM
            CALL SBMMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(27) !*MUSCL_THETA
            CALL THETAMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(28) !*MUSCL_SBM standard
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(2),tau_input=Tau_Standard_Limiting,theta_input=Theta_Standard_Limiting)
          CASE(29) !*MUSCL_SBM overcompressive
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(2),tau_input=Tau_Overcompressive,theta_input=Theta_Overcompressive)
          CASE(3,4,5,7)
            !*NB: Here I'm entering with conserved variables
            CALL WENO_YDIR(&
                      WsupportBox(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts),& !*Entering with Cons
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(2),&
                      Reconstruction,&
                      RRX(1:nVar,1:nVar,ii,jj),&
                      LLX(1:nVar,1:nVar,ii,jj),&
                      RRY(1:nVar,1:nVar,ii,jj),&
                      LLY(1:nVar,1:nVar,ii,jj))
          CASE DEFAULT
            ErrorMessage = "Reconstruction not implemented in ReconstructionY_INPUT"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT

      !*We transform back
        SELECT CASE (ReconstructedVariable)
          CASE(0) !*Conservative

            !*Do nothing

          CASE(1) !*Characteristic

            SELECT CASE (ABS(Reconstruction))
              CASE(1,10,2,20,21,22,23,24,25,26,27,28,29)

                DO iGP=1,nGPs
                  WM(1:nVar,iGP,ii,jj)=MATMUL(RRY(:,:,ii,jj),WM(1:nVar,iGP,ii,jj))
                  WP(1:nVar,iGP,ii,jj)=MATMUL(RRY(:,:,ii,jj),WP(1:nVar,iGP,ii,jj))
                END DO

              CASE(3,4,5,7)

                !*Do nothing

              CASE DEFAULT

                PRINT*, "YOU SHOULD NOT BE HERE IN ReconstructionY_INPUT_CONSERVED"
                STOP

            END SELECT


          CASE(2) !*Primitive

            DO iGP=1,nGPs
              WM(1:nVar,iGP,ii,jj)=MATMUL(MPtoC(:,:,ii,jj),WM(1:nVar,iGP,ii,jj))
              WP(1:nVar,iGP,ii,jj)=MATMUL(MPtoC(:,:,ii,jj),WP(1:nVar,iGP,ii,jj))
            END DO


          CASE DEFAULT
            ErrorMessage = "Error in ReconstructionY_INPUT_CONSERVED. ReconstructedVariable not available"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT





#ifdef PWCINTROUBLEDCELLS
#ifdef MULTIFLUID
        IF(PRESENT(Troubled_Cell_W)) THEN
          IF (Troubled_Cell_W(ii,jj) .NE. 0) THEN !*REPLACE BY PWC
            WM(1:nVar,nGPs,ii,jj) = W(1:nVar,ii,jj)
            WP(1:nVar,nGPs,ii,jj) = W(1:nVar,ii,jj)
          END IF
        END IF
#else
        PRINT*, "Impossible to use PWCINTROUBLEDCELLS"
        PRINT*, "Disable this flag"
        PRINT*, "You should not be here"
        STOP
#endif
#endif

#ifdef MULTIFLUID
        WM(nVar+1,1,ii,jj)=W(nVar,ii,jj)
        WP(nVar+1,1,ii,jj)=W(nVar,ii,jj)
#endif


      END DO
    END DO

  CASE(-1)

    PRINT*, "Reconstruction of Characteristic variables of Alex' and Shaoshuai not set up for subroutine ReconstructionY_INPUT_CONSERVED"
    STOP

  CASE DEFAULT
    ErrorMessage = "Error in ReconstructionY_INPUT. ReconstructedVariable not available"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT




!-------------------------------------------------------------------------------!
END SUBROUTINE ReconstructionY_INPUT_CONSERVED
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ReconstructionFixX_INPUT_PRIMITIVE(W,vbXs,vbXe,vbYs,vbYe,lbXs,lbXe,lbYs,lbYe,Troubled_Cell_W,ReconstructionInput)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructionFix
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructedVariable
USE MOD_FiniteVolume2D_vars,ONLY: LLX
USE MOD_FiniteVolume2D_vars,ONLY: RRX
USE MOD_FiniteVolume2D_vars,ONLY: LLY
USE MOD_FiniteVolume2D_vars,ONLY: RRY
USE MOD_FiniteVolume2D_vars,ONLY: MCtoP
USE MOD_FiniteVolume2D_vars,ONLY: MPtoC
USE MOD_FiniteVolume2D_vars,ONLY: Tau_Standard_Limiting
USE MOD_FiniteVolume2D_vars,ONLY: Theta_Standard_Limiting
USE MOD_FiniteVolume2D_vars,ONLY: Tau_Overcompressive
USE MOD_FiniteVolume2D_vars,ONLY: Theta_Overcompressive
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,vbXs:vbXe,vbYs:vbYe),    INTENT(IN)          :: W
INTEGER,                                        INTENT(IN)          :: vbXs,vbXe,vbYs,vbYe !*Bounds for the vector
INTEGER,                                        INTENT(IN)          :: lbXs,lbXe,lbYs,lbYe !*Bounds for the loop   
INTEGER, DIMENSION(vbXs:vbXe,vbYs:vbYe),        INTENT(IN)          :: Troubled_Cell_W
INTEGER,                                        INTENT(IN),OPTIONAL :: ReconstructionInput
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj, kk, iGP, indi, indj
INTEGER            :: ReconstructionUsed
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: Wsupport(1:nVar,-nGhosts:nGhosts)
REAL               :: WsupportBox(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts)
REAL               :: LL(1:nVar,1:nVar)
REAL               :: RR(1:nVar,1:nVar)
REAL               :: MM(1:nVar,1:nVar)
INTEGER            :: i_interface
#ifdef MULTIFLUID
REAL               ::   roM,    roP,   roHAT
REAL               ::    uM,     uP,    uHAT
REAL               ::  phiM,   phiP,  phiHAT
REAL               :: pinfM,  pinfP, pinfHAT
REAL               ::    pM,     pP,    pHAT
REAL               ::  gmmM,   gmmP,  gmmHAT
REAL               ::    cM,     cP,    cHAT
REAL               :: Wdummy(1:nVar,1:nGPs)
#endif
!-------------------------------------------------------------------------------!

IF (PRESENT(ReconstructionInput)) THEN
  ReconstructionUsed=ReconstructionInput
ELSE
  ReconstructionUsed=ReconstructionFix
END IF

SELECT CASE (ReconstructedVariable)
  CASE(0,1,2)

    DO jj=lbYs,lbYe
      DO ii=lbXs,lbXe

        IF ( Troubled_Cell_W(ii,jj) .EQ. 0) THEN
          CYCLE
        END IF

        SELECT CASE (ReconstructedVariable)
          CASE(0) !*Conservative
            DO kk=-nGhosts,nGhosts
              Wsupport(1:nVar,kk)=MATMUL(MPtoC(:,:,ii,jj),W(1:nVar,ii+kk,jj)) !*Prim->Cons
            END DO
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=MATMUL(MPtoC(:,:,ii,jj),W(1:nVar,ii+indi,jj+indj))  !*Prim->Cons
              END DO
            END DO
          CASE(1) !*Characteristic
            DO kk=-nGhosts,nGhosts
              Wsupport(1:nVar,kk)=MATMUL( LLX(:,:,ii,jj),W(1:nVar,ii+kk,jj) ) !*Prim->Char
            END DO
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=W(1:nVar,ii+indi,jj+indj)  !*Prim
              END DO
            END DO
          CASE(2) !*Primitive
            Wsupport(1:nVar,-nGhosts:nGhosts)=W(1:nVar,-nGhosts+ii:ii+nGhosts,jj)
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=W(1:nVar,ii+indi,jj+indj)
              END DO
            END DO
          CASE DEFAULT
            ErrorMessage = "Error in ReconstructionFixX. ReconstructedVariable not available"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT


        SELECT CASE (ABS(ReconstructionUsed))
          CASE(1,10)
            WM(1:nVar,nGPs,ii,jj) = Wsupport(1:nVar,0)
            WP(1:nVar,nGPs,ii,jj) = Wsupport(1:nVar,0)
          CASE(2,20) !*MUSCL
            CALL MUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(21) !*k-MUSCL
            CALL kMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(22) !*MUSCL_CO
            CALL COMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(23) !*MUSCL_VL
            CALL VLMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(24) !*MUSCL_M
            CALL M2MUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(25) !*MUSCL_VA
            CALL VAMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(26) !*MUSCL_SBM
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(27) !*MUSCL_THETA
            CALL THETAMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(28) !*MUSCL_SBM standard
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1),tau_input=Tau_Standard_Limiting, theta_input=Theta_Standard_Limiting)
          CASE(29) !*MUSCL_SBM overcompressive
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1),tau_input=Tau_Overcompressive, theta_input=Theta_Overcompressive)
          CASE(3,4,5,7)
            !*NB: Here I'm entering with conserved variables
            CALL WENO_XDIR(&
                        WsupportBox(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts),& !*Entering with Cons
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(1),&
                        ReconstructionUsed,&
                        RRX(1:nVar,1:nVar,ii,jj),&
                        LLX(1:nVar,1:nVar,ii,jj),&
                        RRY(1:nVar,1:nVar,ii,jj),&
                        LLY(1:nVar,1:nVar,ii,jj))
          CASE DEFAULT
            ErrorMessage = "ReconstructionFixX_INPUT not implemented"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT


      !*We transform back
        SELECT CASE (ReconstructedVariable)
          CASE(0) !*Conservative

            DO iGP=1,nGPs
              WM(1:nVar,iGP,ii,jj)=MATMUL(MCtoP(:,:,ii,jj),WM(1:nVar,iGP,ii,jj))
              WP(1:nVar,iGP,ii,jj)=MATMUL(MCtoP(:,:,ii,jj),WP(1:nVar,iGP,ii,jj))
            END DO

          CASE(1) !*Characteristic

            SELECT CASE (ABS(ReconstructionUsed))
              CASE(1,10,2,20,21,22,23,24,25,26,27,28,29)

                DO iGP=1,nGPs
                  WM(1:nVar,iGP,ii,jj)=MATMUL(RRX(:,:,ii,jj),WM(1:nVar,iGP,ii,jj))
                  WP(1:nVar,iGP,ii,jj)=MATMUL(RRX(:,:,ii,jj),WP(1:nVar,iGP,ii,jj))
                END DO

              CASE(3,4,5,7)

                !*DO iGP=1,nGPs
                !*  WM(1:nVar,iGP,ii,jj)=WM(1:nVar,iGP,ii,jj)
                !*  WP(1:nVar,iGP,ii,jj)=WP(1:nVar,iGP,ii,jj)
                !*END DO

              CASE DEFAULT

                PRINT*, "YOU SHOULD NOT BE HERE IN ReconstructionFixX_INPUT_PRIMITIVE"
                STOP

            END SELECT


          CASE(2) !*Primitive

            !*Do nothing, primitive variables have been reconstructed

          CASE DEFAULT
            ErrorMessage = "Error in ReconstructionFixX_INPUT_PRIMITIVE. ReconstructedVariable not available"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT




#ifdef MULTIFLUID
        WM(nVar+1,1,ii,jj)=W(nVar,ii,jj)
        WP(nVar+1,1,ii,jj)=W(nVar,ii,jj)
#endif


  END DO
END DO




  CASE(-1)

#ifndef MULTIFLUID
    PRINT*, "Reconstruction of Characteristic variables not set up for MULTIFLUID not active"
    STOP
#else

#if defined(MOMENTUMINPRIMITIVEVARIABLES) && !defined(RECONSTRUCTIONVELOCITY)
    PRINT*, "Reconstruction of Characteristic variables not set up for MOMENTUMINPRIMITIVEVARIABLES without RECONSTRUCTIONVELOCITY"
    STOP
#endif

    DO jj=lbYs,lbYe
      DO ii=lbXs,lbXe

        IF ( Troubled_Cell_W(ii,jj) .EQ. 0) THEN
          CYCLE
        END IF      

        !*LEFT INTERFACE
        i_interface=ii-1

        roM=W(1,i_interface,jj)
        roP=W(1,i_interface+1,jj)  
        !*NB: I always enter with the velocity here
        uM =W(2,i_interface,jj)
        uP =W(2,i_interface+1,jj)  
        !*NB: I always enter with rhophi here
        phiM =W(nVar,i_interface,jj)/W(1,i_interface,jj)
        phiP =W(nVar,i_interface+1,jj)/W(1,i_interface+1,jj)   
        GmmM=Get_Gamma(W(nVar,i_interface,jj))
        pinfM=Get_Pressure_Infinity(W(nVar,i_interface,jj))
        GmmP=Get_Gamma(W(nVar,i_interface+1,jj))
        pinfP=Get_Pressure_Infinity(W(nVar,i_interface+1,jj))
        pM=W(4,i_interface,jj)
        pP=W(4,i_interface+1,jj)  

        !*Roe averages 
        roHAT=SQRT(roM*roP) 
        uHAT =(SQRT(roM)*uM+SQRT(roP)*uP)/(SQRT(roM)+SQRT(roP))                   
        phiHAT=(SQRT(roM)*phiM+SQRT(roP)*phiP)/(SQRT(roM)+SQRT(roP))                    
        pinfHAT=(SQRT(roM)*pinfM+SQRT(roP)*pinfP)/(SQRT(roM)+SQRT(roP))                     
        pHat=(SQRT(roM)*pM+SQRT(roP)*pP)/(SQRT(roM)+SQRT(roP))                     
        gmmHAT=(SQRT(roM)*gmmM+SQRT(roP)*gmmP)/(SQRT(roM)+SQRT(roP))                   
        cHat=SQRT(ABS(gmmHAT*(pHat+pinfHAT)/roHAT)) 

        LL=0.0
        RR=0.0

        RR(1,1)=1.0/cHAT**2
        RR(2,1)=-1.0/(roHat*cHAT)
        RR(4,1)=1.0

        RR(nVar,2)=1.0

        RR(3,3)=1.0

        RR(1,4)=1.0

        RR(1,nVar)=1.0/cHAT**2
        RR(2,nVar)=1.0/(roHat*cHAT)
        RR(4,nVar)=1.0


        LL(4,1)=1.0

        LL(1,2)=-roHAT*cHAT/2.0
        LL(nVar,2)=roHAT*cHAT/2.0

        LL(3,3)=1.0

        LL(1,4)=1.0/2.0
        LL(4,4)=-1.0/cHAT**2
        LL(nVar,4)=1.0/2.0

        LL(2,nVar)=1.0

        DO kk=-nGhosts,nGhosts
          !*NB: I need to have phi in the last component
          Wsupport(1:nVar,kk)=W(1:nVar,ii+kk,jj)
          Wsupport(nVar,kk)=Wsupport(nVar,kk)/Wsupport(1,kk)
          Wsupport(1:nVar,kk)=MATMUL(LL(:,:),Wsupport(1:nVar,kk))
        END DO

        ! MM=MATMUL(LL,RR)
        ! DO indi=1,nVar
        !   PRINT*, indi, MM(indi,:)
        ! END DO
        ! PRINT*


        SELECT CASE (ABS(ReconstructionUsed))
          CASE(1,10)
            WM(1:nVar,nGPs,ii,jj) = Wsupport(1:nVar,0)
          CASE(2,20) !*MUSCL
            CALL MUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      Wdummy(1:nVar,1:nGPs),&
                      MESH_DX(1))
          CASE(21) !*k-MUSCL
            CALL kMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      Wdummy(1:nVar,1:nGPs),&
                      MESH_DX(1))
          CASE(22) !*MUSCL_CO
            CALL COMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      Wdummy(1:nVar,1:nGPs),&
                      MESH_DX(1))
          CASE(23) !*MUSCL_VL
            CALL VLMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      Wdummy(1:nVar,1:nGPs),&
                      MESH_DX(1))
          CASE(24) !*MUSCL_M
            CALL M2MUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      Wdummy(1:nVar,1:nGPs),&
                      MESH_DX(1))
          CASE(25) !*MUSCL_VA
            CALL VAMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      Wdummy(1:nVar,1:nGPs),&
                      MESH_DX(1))
          CASE(26) !*MUSCL_SBM
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      Wdummy(1:nVar,1:nGPs),&
                      MESH_DX(1))
          CASE(27) !*MUSCL_THETA
            CALL THETAMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      Wdummy(1:nVar,1:nGPs),&
                      MESH_DX(1))
          CASE(28) !*MUSCL_SBM standard
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      Wdummy(1:nVar,1:nGPs),&
                      MESH_DX(1),tau_input=Tau_Standard_Limiting,theta_input=Theta_Standard_Limiting)
          CASE(29) !*MUSCL_SBM overcompressive
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      Wdummy(1:nVar,1:nGPs),&
                      MESH_DX(1),tau_input=Tau_Overcompressive,theta_input=Theta_Overcompressive)
          CASE(3,4,5,7)
            !*NB: No shock indicator
            PRINT*, "WENO reconstruction not allowed when reconstructing Chracteristic variables in Active Flux"
            STOP
          CASE DEFAULT
            ErrorMessage = "ReconstructionFixX_INPUT not implemented"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT

        DO iGP=1,nGPs
          WM(1:nVar,iGP,ii,jj)=MATMUL(RR(:,:),WM(1:nVar,iGP,ii,jj))
          WM(nVar,iGP,ii,jj)=WM(nVar,iGP,ii,jj)*WM(1,iGP,ii,jj)
        END DO

        !*RIGHT INTERFACE
        i_interface=ii

        roM=W(1,i_interface,jj)
        roP=W(1,i_interface+1,jj)  

        !*NB: I always enter with the velocity here
        uM =W(2,i_interface,jj)
        uP =W(2,i_interface+1,jj)  

        phiM =W(nVar,i_interface,jj)/W(1,i_interface,jj)
        phiP =W(nVar,i_interface+1,jj)/W(1,i_interface+1,jj)   
        GmmM=Get_Gamma(W(nVar,i_interface,jj))
        pinfM=Get_Pressure_Infinity(W(nVar,i_interface,jj))
        GmmP=Get_Gamma(W(nVar,i_interface+1,jj))
        pinfP=Get_Pressure_Infinity(W(nVar,i_interface+1,jj))
        pM=W(4,i_interface,jj)
        pP=W(4,i_interface+1,jj)  

        !*Roe averages 
        roHAT=SQRT(roM*roP) 
        uHAT =(SQRT(roM)*uM+SQRT(roP)*uP)/(SQRT(roM)+SQRT(roP))                   
        phiHAT=(SQRT(roM)*phiM+SQRT(roP)*phiP)/(SQRT(roM)+SQRT(roP))                    
        pinfHAT=(SQRT(roM)*pinfM+SQRT(roP)*pinfP)/(SQRT(roM)+SQRT(roP))                     
        pHat=(SQRT(roM)*pM+SQRT(roP)*pP)/(SQRT(roM)+SQRT(roP))                     
        gmmHAT=(SQRT(roM)*gmmM+SQRT(roP)*gmmP)/(SQRT(roM)+SQRT(roP))                   
        cHat=SQRT(ABS(gmmHAT*(pHat+pinfHAT)/roHAT)) 

        LL=0.0
        RR=0.0

        RR(1,1)=1.0/cHAT**2
        RR(2,1)=-1.0/(roHat*cHAT)
        RR(4,1)=1.0

        RR(nVar,2)=1.0

        RR(3,3)=1.0

        RR(1,4)=1.0

        RR(1,nVar)=1.0/cHAT**2
        RR(2,nVar)=1.0/(roHat*cHAT)
        RR(4,nVar)=1.0


        LL(4,1)=1.0

        LL(1,2)=-roHAT*cHAT/2.0
        LL(nVar,2)=roHAT*cHAT/2.0

        LL(3,3)=1.0

        LL(1,4)=1.0/2.0
        LL(4,4)=-1.0/cHAT**2
        LL(nVar,4)=1.0/2.0

        LL(2,nVar)=1.0

        DO kk=-nGhosts,nGhosts
          !*NB: I need to have phi in the last component
          Wsupport(1:nVar,kk)=W(1:nVar,ii+kk,jj)
          Wsupport(nVar,kk)=Wsupport(nVar,kk)/Wsupport(1,kk)
          Wsupport(1:nVar,kk)=MATMUL(LL(:,:),Wsupport(1:nVar,kk))
        END DO


        SELECT CASE (ABS(ReconstructionUsed))
          CASE(1,10)
            WP(1:nVar,nGPs,ii,jj) = Wsupport(1:nVar,0)
          CASE(2,20) !*MUSCL
            CALL MUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      Wdummy(1:nVar,1:nGPs),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(21) !*k-MUSCL
            CALL kMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      Wdummy(1:nVar,1:nGPs),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(22) !*MUSCL_CO
            CALL COMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      Wdummy(1:nVar,1:nGPs),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(23) !*MUSCL_VL
            CALL VLMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      Wdummy(1:nVar,1:nGPs),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(24) !*MUSCL_M
            CALL M2MUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      Wdummy(1:nVar,1:nGPs),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(25) !*MUSCL_VA
            CALL VAMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      Wdummy(1:nVar,1:nGPs),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(26) !*MUSCL_SBM
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      Wdummy(1:nVar,1:nGPs),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(27) !*MUSCL_THETA
            CALL THETAMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      Wdummy(1:nVar,1:nGPs),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(28) !*MUSCL_SBM standard
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      Wdummy(1:nVar,1:nGPs),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1),tau_input=Tau_Standard_Limiting,theta_input=Theta_Standard_Limiting)
          CASE(29) !*MUSCL_SBM overcompressive
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      Wdummy(1:nVar,1:nGPs),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1),tau_input=Tau_Overcompressive,theta_input=Theta_Overcompressive)
          CASE(3,4,5,7)
            !*NB: No shock indicator
            PRINT*, "WENO reconstruction not allowed when reconstructing Chracteristic variables in Active Flux"
            STOP
          CASE DEFAULT
            ErrorMessage = "ReconstructionFixX_INPUT not implemented"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT


        DO iGP=1,nGPs
          WP(1:nVar,iGP,ii,jj)=MATMUL(RR(:,:),WP(1:nVar,iGP,ii,jj))
          WP(nVar,iGP,ii,jj)=WP(nVar,iGP,ii,jj)*WP(1,iGP,ii,jj)
        END DO






#ifdef MULTIFLUID
        WM(nVar+1,1,ii,jj)=W(nVar,ii,jj)
        WP(nVar+1,1,ii,jj)=W(nVar,ii,jj)
#endif


  END DO
END DO


#endif

  CASE DEFAULT
    ErrorMessage = "Error in ReconstructionFixX_INPUT. ReconstructedVariable not available"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT




!-------------------------------------------------------------------------------!
END SUBROUTINE ReconstructionFixX_INPUT_PRIMITIVE
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ReconstructionFixY_INPUT_PRIMITIVE(W,vbXs,vbXe,vbYs,vbYe,lbXs,lbXe,lbYs,lbYe,Troubled_Cell_W,ReconstructionInput)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructionFix
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructedVariable
USE MOD_FiniteVolume2D_vars,ONLY: LLX
USE MOD_FiniteVolume2D_vars,ONLY: RRX
USE MOD_FiniteVolume2D_vars,ONLY: LLY
USE MOD_FiniteVolume2D_vars,ONLY: RRY
USE MOD_FiniteVolume2D_vars,ONLY: MCtoP
USE MOD_FiniteVolume2D_vars,ONLY: MPtoC
USE MOD_FiniteVolume2D_vars,ONLY: Tau_Standard_Limiting
USE MOD_FiniteVolume2D_vars,ONLY: Theta_Standard_Limiting
USE MOD_FiniteVolume2D_vars,ONLY: Tau_Overcompressive
USE MOD_FiniteVolume2D_vars,ONLY: Theta_Overcompressive
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,vbXs:vbXe,vbYs:vbYe),    INTENT(IN)          :: W
INTEGER,                                        INTENT(IN)          :: vbXs,vbXe,vbYs,vbYe !*Bounds for the vector
INTEGER,                                        INTENT(IN)          :: lbXs,lbXe,lbYs,lbYe !*Bounds for the loop    
INTEGER, DIMENSION(vbXs:vbXe,vbYs:vbYe),        INTENT(IN)          :: Troubled_Cell_W                        
INTEGER,                                        INTENT(IN),OPTIONAL :: ReconstructionInput
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj, kk, iGP
INTEGER            :: ReconstructionUsed
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: Wsupport(1:nVar,-nGhosts:nGhosts)
REAL               :: WsupportBox(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts)
REAL               :: LL(1:nVar,1:nVar)
REAL               :: RR(1:nVar,1:nVar)
REAL               :: MM(1:nVar,1:nVar)
INTEGER            :: j_interface
REAL               :: a_support, b_support
#ifdef MULTIFLUID
REAL               ::   roM,    roP,   roHAT
REAL               ::    uM,     uP,    uHAT
REAL               ::  phiM,   phiP,  phiHAT
REAL               :: pinfM,  pinfP, pinfHAT
REAL               ::    pM,     pP,    pHAT
REAL               ::  gmmM,   gmmP,  gmmHAT
REAL               ::    cM,     cP,    cHAT
REAL               :: Wdummy(1:nVar,1:nGPs)
#endif
INTEGER            :: indi, indj
!-------------------------------------------------------------------------------!

IF (PRESENT(ReconstructionInput)) THEN
  ReconstructionUsed=ReconstructionInput
ELSE
  ReconstructionUsed=ReconstructionFix
END IF


SELECT CASE (ReconstructedVariable)
  CASE(0,1,2)

    DO jj=lbYs,lbYe
      DO ii=lbXs,lbXe

        IF ( Troubled_Cell_W(ii,jj) .EQ. 0) THEN
          CYCLE
        END IF

        SELECT CASE (ReconstructedVariable)
          CASE(0) !*Conservative
            DO kk=-nGhosts,nGhosts
              Wsupport(1:nVar,kk)=MATMUL(MPtoC(:,:,ii,jj),W(1:nVar,ii,jj+kk)) !*Prim->Cons
            END DO
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=MATMUL(MPtoC(:,:,ii,jj),W(1:nVar,ii+indi,jj+indj))  !*Prim->Cons
              END DO
            END DO
          CASE(1) !*Characteristic
            DO kk=-nGhosts,nGhosts
              Wsupport(1:nVar,kk)=MATMUL( LLY(:,:,ii,jj),W(1:nVar,ii,jj+kk))  !*Prim->Char
            END DO
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=W(1:nVar,ii+indi,jj+indj)  !*Prim->Cons
              END DO
            END DO
          CASE(2) !*Primitive
            Wsupport(1:nVar,-nGhosts:nGhosts)=W(1:nVar,ii,-nGhosts+jj:jj+nGhosts)
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=W(1:nVar,ii+indi,jj+indj)
              END DO
            END DO
          CASE DEFAULT
            ErrorMessage = "Error in ReconstructionFixY. ReconstructedVariable not available"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT



        SELECT CASE (ABS(ReconstructionUsed))
          CASE(1,10)
            WM(1:nVar,nGPs,ii,jj) = Wsupport(1:nVar,0)
            WP(1:nVar,nGPs,ii,jj) = Wsupport(1:nVar,0)
          CASE(2,20) !*MUSCL
            CALL MUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(21) !*k-MUSCL
            CALL kMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(22) !*MUSCL_CO
            CALL COMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(23) !*MUSCL_VL
            CALL VLMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(24) !*MUSCL_M
            CALL M2MUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(25) !*MUSCL_VA
            CALL VAMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(26) !*MUSCL_SBM
            CALL SBMMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(27) !*MUSCL_THETA
            CALL THETAMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(28) !*MUSCL_SBM standard
            CALL SBMMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2),tau_input=Tau_Standard_Limiting,theta_input=Theta_Standard_Limiting)
          CASE(29) !*MUSCL_SBM overcompressive
            CALL SBMMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2),tau_input=Tau_Overcompressive,theta_input=Theta_Overcompressive)
          CASE(3,4,5,7)
            !*NB: Here I'm entering with conserved variables
            CALL WENO_YDIR(&
                      WsupportBox(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts),& !*Entering with Cons
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(2),&
                      ReconstructionUsed,&
                      RRX(1:nVar,1:nVar,ii,jj),&
                      LLX(1:nVar,1:nVar,ii,jj),&
                      RRY(1:nVar,1:nVar,ii,jj),&
                      LLY(1:nVar,1:nVar,ii,jj))
          CASE DEFAULT
            ErrorMessage = "Reconstruction not implemented in ReconstructionFixY_INPUT"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT

      !*We transform back
        SELECT CASE (ReconstructedVariable)
          CASE(0) !*Conservative

            DO iGP=1,nGPs
              WM(1:nVar,iGP,ii,jj)=MATMUL(MCtoP(:,:,ii,jj),WM(1:nVar,iGP,ii,jj))
              WP(1:nVar,iGP,ii,jj)=MATMUL(MCtoP(:,:,ii,jj),WP(1:nVar,iGP,ii,jj))
            END DO

          CASE(1) !*Characteristic

            SELECT CASE (ABS(ReconstructionUsed))
              CASE(1,10,2,20,21,22,23,24,25,26,27,28,29)

                DO iGP=1,nGPs
                  WM(1:nVar,iGP,ii,jj)=MATMUL(RRY(:,:,ii,jj),WM(1:nVar,iGP,ii,jj))
                  WP(1:nVar,iGP,ii,jj)=MATMUL(RRY(:,:,ii,jj),WP(1:nVar,iGP,ii,jj))
                END DO

              CASE(3,4,5,7)

                !*DO iGP=1,nGPs
                !*  WM(1:nVar,iGP,ii,jj)=WM(1:nVar,iGP,ii,jj)
                !*  WP(1:nVar,iGP,ii,jj)=WP(1:nVar,iGP,ii,jj)
                !*END DO

              CASE DEFAULT

                PRINT*, "YOU SHOULD NOT BE HERE IN ReconstructionFixY_INPUT_PRIMITIVE"
                STOP

            END SELECT


          CASE(2) !*Primitive

            !*Do nothing, primitive variables have been reconstructed

          CASE DEFAULT
            ErrorMessage = "Error in ReconstructionFixY_INPUT_PRIMITIVE. ReconstructedVariable not available"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT





#ifdef MULTIFLUID
        WM(nVar+1,1,ii,jj)=W(nVar,ii,jj)
        WP(nVar+1,1,ii,jj)=W(nVar,ii,jj)
#endif


      END DO
    END DO

  CASE(-1)

#ifndef MULTIFLUID
    PRINT*, "Reconstruction of Characteristic variables not set up for MULTIFLUID not active"
    STOP
#else

#if defined(MOMENTUMINPRIMITIVEVARIABLES) && !defined(RECONSTRUCTIONVELOCITY)
    PRINT*, "Reconstruction of Characteristic variables not set up for MOMENTUMINPRIMITIVEVARIABLES without RECONSTRUCTIONVELOCITY"
    STOP
#endif

    DO jj=lbYs,lbYe
      DO ii=lbXs,lbXe

        IF ( Troubled_Cell_W(ii,jj) .EQ. 0) THEN
          CYCLE
        END IF

        !*BOTTOM INTERFACE
        j_interface=jj-1

        roM=W(1,ii,j_interface)
        roP=W(1,ii,j_interface+1)  
        !*NB: I always enter with the velocity here
        uM =W(3,ii,j_interface)   !*NB: This time I take the v component
        uP =W(3,ii,j_interface+1) !*NB: This time I take the v component 
        !*NB: I always enter with rhophi here
        phiM =W(nVar,ii,j_interface)/W(1,ii,j_interface)
        phiP =W(nVar,ii,j_interface+1)/W(1,ii,j_interface+1)   
        GmmM=Get_Gamma(W(nVar,ii,j_interface))
        pinfM=Get_Pressure_Infinity(W(nVar,ii,j_interface))
        GmmP=Get_Gamma(W(nVar,ii,j_interface+1))
        pinfP=Get_Pressure_Infinity(W(nVar,ii,j_interface+1))
        pM=W(4,ii,j_interface)
        pP=W(4,ii,j_interface+1)  

        !*Roe averages 
        roHAT=SQRT(roM*roP) 
        uHAT =(SQRT(roM)*uM+SQRT(roP)*uP)/(SQRT(roM)+SQRT(roP))                   
        phiHAT=(SQRT(roM)*phiM+SQRT(roP)*phiP)/(SQRT(roM)+SQRT(roP))                    
        pinfHAT=(SQRT(roM)*pinfM+SQRT(roP)*pinfP)/(SQRT(roM)+SQRT(roP))                     
        pHat=(SQRT(roM)*pM+SQRT(roP)*pP)/(SQRT(roM)+SQRT(roP))                     
        gmmHAT=(SQRT(roM)*gmmM+SQRT(roP)*gmmP)/(SQRT(roM)+SQRT(roP))                   
        cHat=SQRT(ABS(gmmHAT*(pHat+pinfHAT)/roHAT)) 

        LL=0.0
        RR=0.0

        RR(1,1)=1.0/cHAT**2
        RR(2,1)=-1.0/(roHat*cHAT)
        RR(4,1)=1.0

        RR(nVar,2)=1.0

        RR(3,3)=1.0

        RR(1,4)=1.0

        RR(1,nVar)=1.0/cHAT**2
        RR(2,nVar)=1.0/(roHat*cHAT)
        RR(4,nVar)=1.0


        LL(4,1)=1.0

        LL(1,2)=-roHAT*cHAT/2.0
        LL(nVar,2)=roHAT*cHAT/2.0

        LL(3,3)=1.0

        LL(1,4)=1.0/2.0
        LL(4,4)=-1.0/cHAT**2
        LL(nVar,4)=1.0/2.0

        LL(2,nVar)=1.0

        DO kk=-nGhosts,nGhosts
          Wsupport(1:nVar,kk)=W(1:nVar,ii,jj+kk)
          !*NB: Switching velocities
          Wsupport(2,kk)=W(3,ii,jj+kk) 
          Wsupport(3,kk)=W(2,ii,jj+kk) 
          !*NB: I need to have phi in the last component
          Wsupport(nVar,kk)=Wsupport(nVar,kk)/Wsupport(1,kk)
          Wsupport(1:nVar,kk)=MATMUL(LL(:,:),Wsupport(1:nVar,kk))
        END DO



        SELECT CASE (ABS(ReconstructionUsed))
          CASE(1,10)
            WM(1:nVar,nGPs,ii,jj) = Wsupport(1:nVar,0)
          CASE(2,20) !*MUSCL
            CALL MUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        Wdummy(1:nVar,1:nGPs),&
                        MESH_DX(2))
          CASE(21) !*k-MUSCL
            CALL kMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        Wdummy(1:nVar,1:nGPs),&
                        MESH_DX(2))
          CASE(22) !*MUSCL_CO
            CALL COMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        Wdummy(1:nVar,1:nGPs),&
                        MESH_DX(2))
          CASE(23) !*MUSCL_VL
            CALL VLMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        Wdummy(1:nVar,1:nGPs),&
                        MESH_DX(2))
          CASE(24) !*MUSCL_M
            CALL M2MUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        Wdummy(1:nVar,1:nGPs),&
                        MESH_DX(2))
          CASE(25) !*MUSCL_VA
            CALL VAMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        Wdummy(1:nVar,1:nGPs),&
                        MESH_DX(2))
          CASE(26) !*MUSCL_SBM
            CALL SBMMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        Wdummy(1:nVar,1:nGPs),&
                        MESH_DX(2))
          CASE(27) !*MUSCL_THETA
            CALL THETAMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        Wdummy(1:nVar,1:nGPs),&
                        MESH_DX(2))
          CASE(28) !*MUSCL_SBM standard
            CALL SBMMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        Wdummy(1:nVar,1:nGPs),&
                        MESH_DX(2),tau_input=Tau_Standard_Limiting,theta_input=Theta_Standard_Limiting)
          CASE(29) !*MUSCL_SBM overcompressive
            CALL SBMMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        Wdummy(1:nVar,1:nGPs),&
                        MESH_DX(2),tau_input=Tau_Overcompressive,theta_input=Theta_Overcompressive)
          CASE(3,4,5,7)
            !*NB: No shock indicator
            PRINT*, "WENO reconstruction not allowed when reconstructing Chracteristic variables in Active Flux"
            STOP
          CASE DEFAULT
            ErrorMessage = "Reconstruction not implemented in ReconstructionFixY_INPUT"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT

        DO iGP=1,nGPs
          WM(1:nVar,iGP,ii,jj)=MATMUL(RR(:,:),WM(1:nVar,iGP,ii,jj))
          WM(nVar,iGP,ii,jj)=WM(nVar,iGP,ii,jj)*WM(1,iGP,ii,jj)
          !*NB: Switching back the velocities
          a_support=WM(2,iGP,ii,jj) 
          b_support=WM(3,iGP,ii,jj) 
          WM(2,iGP,ii,jj)=b_support 
          WM(3,iGP,ii,jj)=a_support 
        END DO

        !*UPPER INTERFACE
        j_interface=jj

        roM=W(1,ii,j_interface)
        roP=W(1,ii,j_interface+1)  
        !*NB: I always enter with the velocity here
        uM =W(3,ii,j_interface)   !*NB: This time I take the v component
        uP =W(3,ii,j_interface+1) !*NB: This time I take the v component 
        !*NB: I always enter with rhophi here
        phiM =W(nVar,ii,j_interface)/W(1,ii,j_interface)
        phiP =W(nVar,ii,j_interface+1)/W(1,ii,j_interface+1)   
        GmmM=Get_Gamma(W(nVar,ii,j_interface))
        pinfM=Get_Pressure_Infinity(W(nVar,ii,j_interface))
        GmmP=Get_Gamma(W(nVar,ii,j_interface+1))
        pinfP=Get_Pressure_Infinity(W(nVar,ii,j_interface+1))
        pM=W(4,ii,j_interface)
        pP=W(4,ii,j_interface+1)  

        !*Roe averages 
        roHAT=SQRT(roM*roP) 
        uHAT =(SQRT(roM)*uM+SQRT(roP)*uP)/(SQRT(roM)+SQRT(roP))                   
        phiHAT=(SQRT(roM)*phiM+SQRT(roP)*phiP)/(SQRT(roM)+SQRT(roP))                    
        pinfHAT=(SQRT(roM)*pinfM+SQRT(roP)*pinfP)/(SQRT(roM)+SQRT(roP))                     
        pHat=(SQRT(roM)*pM+SQRT(roP)*pP)/(SQRT(roM)+SQRT(roP))                     
        gmmHAT=(SQRT(roM)*gmmM+SQRT(roP)*gmmP)/(SQRT(roM)+SQRT(roP))                   
        cHat=SQRT(ABS(gmmHAT*(pHat+pinfHAT)/roHAT)) 

        LL=0.0
        RR=0.0

        RR(1,1)=1.0/cHAT**2
        RR(2,1)=-1.0/(roHat*cHAT)
        RR(4,1)=1.0

        RR(nVar,2)=1.0

        RR(3,3)=1.0

        RR(1,4)=1.0

        RR(1,nVar)=1.0/cHAT**2
        RR(2,nVar)=1.0/(roHat*cHAT)
        RR(4,nVar)=1.0


        LL(4,1)=1.0

        LL(1,2)=-roHAT*cHAT/2.0
        LL(nVar,2)=roHAT*cHAT/2.0

        LL(3,3)=1.0

        LL(1,4)=1.0/2.0
        LL(4,4)=-1.0/cHAT**2
        LL(nVar,4)=1.0/2.0

        LL(2,nVar)=1.0

        DO kk=-nGhosts,nGhosts
          Wsupport(1:nVar,kk)=W(1:nVar,ii,jj+kk)
          !*NB: Switching velocities
          Wsupport(2,kk)=W(3,ii,jj+kk) 
          Wsupport(3,kk)=W(2,ii,jj+kk) 
          !*NB: I need to have phi in the last component
          Wsupport(nVar,kk)=Wsupport(nVar,kk)/Wsupport(1,kk)
          Wsupport(1:nVar,kk)=MATMUL(LL(:,:),Wsupport(1:nVar,kk))
        END DO

        SELECT CASE (ABS(ReconstructionUsed))
          CASE(1,10)
            WP(1:nVar,nGPs,ii,jj) = Wsupport(1:nVar,0)
          CASE(2,20) !*MUSCL
            CALL MUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        Wdummy(1:nVar,1:nGPs),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(21) !*k-MUSCL
            CALL kMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        Wdummy(1:nVar,1:nGPs),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(22) !*MUSCL_CO
            CALL COMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        Wdummy(1:nVar,1:nGPs),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(23) !*MUSCL_VL
            CALL VLMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        Wdummy(1:nVar,1:nGPs),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(24) !*MUSCL_M
            CALL M2MUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        Wdummy(1:nVar,1:nGPs),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(25) !*MUSCL_VA
            CALL VAMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        Wdummy(1:nVar,1:nGPs),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(26) !*MUSCL_SBM
            CALL SBMMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        Wdummy(1:nVar,1:nGPs),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(27) !*MUSCL_THETA
            CALL THETAMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        Wdummy(1:nVar,1:nGPs),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(28) !*MUSCL_SBM standard
            CALL SBMMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        Wdummy(1:nVar,1:nGPs),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2),tau_input=Tau_Standard_Limiting,theta_input=Theta_Standard_Limiting)
          CASE(29) !*MUSCL_SBM overcompressive
            CALL SBMMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        Wdummy(1:nVar,1:nGPs),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2),tau_input=Tau_Overcompressive,theta_input=Theta_Overcompressive)
          CASE(3,4,5,7)
            !*NB: No shock indicator
            PRINT*, "WENO reconstruction not allowed when reconstructing Chracteristic variables in Active Flux"
            STOP
          CASE DEFAULT
            ErrorMessage = "Reconstruction not implemented in ReconstructionFixY_INPUT"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT



        DO iGP=1,nGPs
          WP(1:nVar,iGP,ii,jj)=MATMUL(RR(:,:),WP(1:nVar,iGP,ii,jj))
          WP(nVar,iGP,ii,jj)=WP(nVar,iGP,ii,jj)*WP(1,iGP,ii,jj)
          !*NB: Switching back the velocities
          a_support=WP(2,iGP,ii,jj) 
          b_support=WP(3,iGP,ii,jj) 
          WP(2,iGP,ii,jj)=b_support 
          WP(3,iGP,ii,jj)=a_support 
        END DO





#ifdef MULTIFLUID
        WM(nVar+1,1,ii,jj)=W(nVar,ii,jj)
        WP(nVar+1,1,ii,jj)=W(nVar,ii,jj)
#endif




  END DO
END DO


#endif

  CASE DEFAULT
    ErrorMessage = "Error in ReconstructionFixY_INPUT. ReconstructedVariable not available"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT




!-------------------------------------------------------------------------------!
END SUBROUTINE ReconstructionFixY_INPUT_PRIMITIVE
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ReconstructionFixX_INPUT_CONSERVED(W,vbXs,vbXe,vbYs,vbYe,lbXs,lbXe,lbYs,lbYe,Troubled_Cell_W)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructionFix
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructedVariable
USE MOD_FiniteVolume2D_vars,ONLY: LLX
USE MOD_FiniteVolume2D_vars,ONLY: RRX
USE MOD_FiniteVolume2D_vars,ONLY: LLY
USE MOD_FiniteVolume2D_vars,ONLY: RRY
USE MOD_FiniteVolume2D_vars,ONLY: MCtoP
USE MOD_FiniteVolume2D_vars,ONLY: MPtoC
USE MOD_FiniteVolume2D_vars,ONLY: Tau_Standard_Limiting
USE MOD_FiniteVolume2D_vars,ONLY: Theta_Standard_Limiting
USE MOD_FiniteVolume2D_vars,ONLY: Tau_Overcompressive
USE MOD_FiniteVolume2D_vars,ONLY: Theta_Overcompressive
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,vbXs:vbXe,vbYs:vbYe),    INTENT(IN)          :: W
INTEGER,                                        INTENT(IN)          :: vbXs,vbXe,vbYs,vbYe !*Bounds for the vector
INTEGER,                                        INTENT(IN)          :: lbXs,lbXe,lbYs,lbYe !*Bounds for the loop   
INTEGER, DIMENSION(vbXs:vbXe,vbYs:vbYe),        INTENT(IN)          :: Troubled_Cell_W                        
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj, kk, iGP, indi, indj
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: Wsupport(1:nVar,-nGhosts:nGhosts)
REAL               :: WsupportBox(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts)
REAL               :: LL(1:nVar,1:nVar)
REAL               :: RR(1:nVar,1:nVar)
REAL               :: MM(1:nVar,1:nVar)
INTEGER            :: i_interface
#ifdef MULTIFLUID
REAL               ::   roM,    roP,   roHAT
REAL               ::    uM,     uP,    uHAT
REAL               ::  phiM,   phiP,  phiHAT
REAL               :: pinfM,  pinfP, pinfHAT
REAL               ::    pM,     pP,    pHAT
REAL               ::  gmmM,   gmmP,  gmmHAT
REAL               ::    cM,     cP,    cHAT
REAL               :: Wdummy(1:nVar,1:nGPs)
#endif
!-------------------------------------------------------------------------------!


SELECT CASE (ReconstructedVariable)
  CASE(0,1,2)

    DO jj=lbYs,lbYe
      DO ii=lbXs,lbXe

        IF ( Troubled_Cell_W(ii,jj) .EQ. 0) THEN
          CYCLE
        END IF

        SELECT CASE (ReconstructedVariable)
          CASE(0) !*Conserved
            Wsupport(1:nVar,-nGhosts:nGhosts)=W(1:nVar,-nGhosts+ii:ii+nGhosts,jj)
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=W(1:nVar,ii+indi,jj+indj)
              END DO
            END DO
          CASE(1) !*Characteristic
            DO kk=-nGhosts,nGhosts
              Wsupport(1:nVar,kk)=MATMUL( LLX(:,:,ii,jj),W(1:nVar,ii+kk,jj) ) !*Cons->Char
            END DO
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=W(1:nVar,ii+indi,jj+indj)  !*Cons
              END DO
            END DO
          CASE(2) !*Primitive
            DO kk=-nGhosts,nGhosts
              Wsupport(1:nVar,kk)=MATMUL(MCtoP(:,:,ii,jj),W(1:nVar,ii+kk,jj)) !*Cons->Prim
            END DO
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=MATMUL(MCtoP(:,:,ii,jj),W(1:nVar,ii+indi,jj+indj))  !*Cons->Prim
              END DO
            END DO
          CASE DEFAULT
            ErrorMessage = "Error in ReconstructionFixX. ReconstructedVariable not available"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT


        SELECT CASE (ABS(ReconstructionFix))
          CASE(1,10)
            WM(1:nVar,nGPs,ii,jj) = Wsupport(1:nVar,0)
            WP(1:nVar,nGPs,ii,jj) = Wsupport(1:nVar,0)
          CASE(2,20) !*MUSCL
            CALL MUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(21) !*k-MUSCL
            CALL kMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(22) !*MUSCL_CO
            CALL COMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(23) !*MUSCL_VL
            CALL VLMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(24) !*MUSCL_M
            CALL M2MUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(25) !*MUSCL_VA
            CALL VAMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(26) !*MUSCL_SBM
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(27) !*MUSCL_THETA
            CALL THETAMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1))
          CASE(28) !*MUSCL_SBM standard
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1),tau_input=Tau_Standard_Limiting,theta_input=Theta_Standard_Limiting)
          CASE(29) !*MUSCL_SBM overcompressive
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(1),tau_input=Tau_Overcompressive,theta_input=Theta_Overcompressive)
          CASE(3,4,5,7)
            !*NB: Here I'm entering with conserved variables
            CALL WENO_XDIR(&
                        WsupportBox(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts),& !*Entering with Cons
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(1),&
                        ReconstructionFix,&
                        RRX(1:nVar,1:nVar,ii,jj),&
                        LLX(1:nVar,1:nVar,ii,jj),&
                        RRY(1:nVar,1:nVar,ii,jj),&
                        LLY(1:nVar,1:nVar,ii,jj))
          CASE DEFAULT
            ErrorMessage = "ReconstructionFixX_INPUT not implemented"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT


      !*We transform back
        SELECT CASE (ReconstructedVariable)
          CASE(0) !*Conservative

            !*Do nothing

          CASE(1) !*Characteristic

            SELECT CASE (ABS(ReconstructionFix))
              CASE(1,10,2,20,21,22,23,24,25,26,27,28,29)

                DO iGP=1,nGPs
                  WM(1:nVar,iGP,ii,jj)=MATMUL(RRX(:,:,ii,jj),WM(1:nVar,iGP,ii,jj))
                  WP(1:nVar,iGP,ii,jj)=MATMUL(RRX(:,:,ii,jj),WP(1:nVar,iGP,ii,jj))
                END DO

              CASE(3,4,5,7)

                !*Do nothing

              CASE DEFAULT

                PRINT*, "YOU SHOULD NOT BE HERE IN ReconstructionFixX_INPUT_CONSERVED"
                STOP

            END SELECT


          CASE(2) !*Primitive

                DO iGP=1,nGPs
                  WM(1:nVar,iGP,ii,jj)=MATMUL(MPtoC(:,:,ii,jj),WM(1:nVar,iGP,ii,jj))
                  WP(1:nVar,iGP,ii,jj)=MATMUL(MPtoC(:,:,ii,jj),WP(1:nVar,iGP,ii,jj))
                END DO

          CASE DEFAULT
            ErrorMessage = "Error in ReconstructionFixX_INPUT_CONSERVED. ReconstructedVariable not available"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT


#ifdef MULTIFLUID
        WM(nVar+1,1,ii,jj)=W(nVar,ii,jj)
        WP(nVar+1,1,ii,jj)=W(nVar,ii,jj)
#endif


  END DO
END DO




  CASE(-1)

    PRINT*, "Reconstruction of Characteristic variables of Alex' and Shaoshuai not set up for subroutine ReconstructionFixX_INPUT_CONSERVED"
    STOP

  CASE DEFAULT
    ErrorMessage = "Error in ReconstructionFixX_INPUT. ReconstructedVariable not available"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT




!-------------------------------------------------------------------------------!
END SUBROUTINE ReconstructionFixX_INPUT_CONSERVED
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ReconstructionFixY_INPUT_CONSERVED(W,vbXs,vbXe,vbYs,vbYe,lbXs,lbXe,lbYs,lbYe,Troubled_Cell_W)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: WM
USE MOD_FiniteVolume2D_vars,ONLY: WP
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructionFix
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructedVariable
USE MOD_FiniteVolume2D_vars,ONLY: LLX
USE MOD_FiniteVolume2D_vars,ONLY: RRX
USE MOD_FiniteVolume2D_vars,ONLY: LLY
USE MOD_FiniteVolume2D_vars,ONLY: RRY
USE MOD_FiniteVolume2D_vars,ONLY: MCtoP
USE MOD_FiniteVolume2D_vars,ONLY: MPtoC
USE MOD_FiniteVolume2D_vars,ONLY: Tau_Standard_Limiting
USE MOD_FiniteVolume2D_vars,ONLY: Theta_Standard_Limiting
USE MOD_FiniteVolume2D_vars,ONLY: Tau_Overcompressive
USE MOD_FiniteVolume2D_vars,ONLY: Theta_Overcompressive
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL, DIMENSION(1:nVar,vbXs:vbXe,vbYs:vbYe),    INTENT(IN)          :: W
INTEGER,                                        INTENT(IN)          :: vbXs,vbXe,vbYs,vbYe !*Bounds for the vector
INTEGER,                                        INTENT(IN)          :: lbXs,lbXe,lbYs,lbYe !*Bounds for the loop    
INTEGER, DIMENSION(vbXs:vbXe,vbYs:vbYe),        INTENT(IN)          :: Troubled_Cell_W                        

!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj, kk, iGP
CHARACTER(LEN=255) :: ErrorMessage
REAL               :: Wsupport(1:nVar,-nGhosts:nGhosts)
REAL               :: WsupportBox(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts)
REAL               :: LL(1:nVar,1:nVar)
REAL               :: RR(1:nVar,1:nVar)
REAL               :: MM(1:nVar,1:nVar)
INTEGER            :: j_interface
REAL               :: a_support, b_support
#ifdef MULTIFLUID
REAL               ::   roM,    roP,   roHAT
REAL               ::    uM,     uP,    uHAT
REAL               ::  phiM,   phiP,  phiHAT
REAL               :: pinfM,  pinfP, pinfHAT
REAL               ::    pM,     pP,    pHAT
REAL               ::  gmmM,   gmmP,  gmmHAT
REAL               ::    cM,     cP,    cHAT
REAL               :: Wdummy(1:nVar,1:nGPs)
#endif
INTEGER            :: indi, indj
!-------------------------------------------------------------------------------!

SELECT CASE (ReconstructedVariable)
  CASE(0,1,2)

    DO jj=lbYs,lbYe
      DO ii=lbXs,lbXe

        IF ( Troubled_Cell_W(ii,jj) .EQ. 0) THEN
          CYCLE
        END IF

        SELECT CASE (ReconstructedVariable)
          CASE(0) !*Conserved
            Wsupport(1:nVar,-nGhosts:nGhosts)=W(1:nVar,ii,-nGhosts+jj:jj+nGhosts)
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=W(1:nVar,ii+indi,jj+indj)
              END DO
            END DO

          CASE(1) !*Characteristic
            DO kk=-nGhosts,nGhosts
              Wsupport(1:nVar,kk)=MATMUL( LLY(:,:,ii,jj),W(1:nVar,ii,jj+kk)) !*Cons->Char
            END DO
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=W(1:nVar,ii+indi,jj+indj) !*Cons
              END DO
            END DO

          CASE(2) !*Primitive
            DO kk=-nGhosts,nGhosts
              Wsupport(1:nVar,kk)=MATMUL(MCtoP(:,:,ii,jj),W(1:nVar,ii,jj+kk)) !*Cons->Prim
            END DO
            DO indj=-nGhosts,nGhosts
              DO indi=-nGhosts,nGhosts
                WsupportBox(1:nVar,indi,indj)=MATMUL(MCtoP(:,:,ii,jj),W(1:nVar,ii+indi,jj+indj))  !*Cons->Prim
              END DO
            END DO


          CASE DEFAULT
            ErrorMessage = "Error in ReconstructionFixY. ReconstructedVariable not available"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT



        SELECT CASE (ABS(ReconstructionFix))
          CASE(1,10)
            WM(1:nVar,nGPs,ii,jj) = Wsupport(1:nVar,0)
            WP(1:nVar,nGPs,ii,jj) = Wsupport(1:nVar,0)
          CASE(2,20) !*MUSCL
            CALL MUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(21) !*k-MUSCL
            CALL kMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(22) !*MUSCL_CO
            CALL COMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(23) !*MUSCL_VL
            CALL VLMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(24) !*MUSCL_M
            CALL M2MUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(25) !*MUSCL_VA
            CALL VAMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(26) !*MUSCL_SBM
            CALL SBMMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(27) !*MUSCL_THETA
            CALL THETAMUSCL(&
                        Wsupport(1:nVar,-nGhosts:nGhosts),&
                        WM(1:nVar,1:nGPs,ii,jj),&
                        WP(1:nVar,1:nGPs,ii,jj),&
                        MESH_DX(2))
          CASE(28) !*MUSCL_SBM standard
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(2),tau_input=Tau_Standard_Limiting,theta_input=Theta_Standard_Limiting)
          CASE(29) !*MUSCL_SBM overcompressive
            CALL SBMMUSCL(&
                      Wsupport(1:nVar,-nGhosts:nGhosts),&
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(2),tau_input=Tau_Overcompressive,theta_input=Theta_Overcompressive)
          CASE(3,4,5,7)
            !*NB: Here I'm entering with conserved variables
            CALL WENO_YDIR(&
                      WsupportBox(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts),& !*Entering with Cons
                      WM(1:nVar,1:nGPs,ii,jj),&
                      WP(1:nVar,1:nGPs,ii,jj),&
                      MESH_DX(2),&
                      ReconstructionFix,&
                      RRX(1:nVar,1:nVar,ii,jj),&
                      LLX(1:nVar,1:nVar,ii,jj),&
                      RRY(1:nVar,1:nVar,ii,jj),&
                      LLY(1:nVar,1:nVar,ii,jj))
          CASE DEFAULT
            ErrorMessage = "Reconstruction not implemented in ReconstructionFixY_INPUT"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT

      !*We transform back
        SELECT CASE (ReconstructedVariable)
          CASE(0) !*Conservative

            !*Do nothing

          CASE(1) !*Characteristic

            SELECT CASE (ABS(ReconstructionFix))
              CASE(1,10,2,20,21,22,23,24,25,26,27,28,29)

                DO iGP=1,nGPs
                  WM(1:nVar,iGP,ii,jj)=MATMUL(RRY(:,:,ii,jj),WM(1:nVar,iGP,ii,jj))
                  WP(1:nVar,iGP,ii,jj)=MATMUL(RRY(:,:,ii,jj),WP(1:nVar,iGP,ii,jj))
                END DO

              CASE(3,4,5,7)

                !*Do nothing

              CASE DEFAULT

                PRINT*, "YOU SHOULD NOT BE HERE IN ReconstructionFixY_INPUT_CONSERVED"
                STOP

            END SELECT


          CASE(2) !*Primitive

            DO iGP=1,nGPs
              WM(1:nVar,iGP,ii,jj)=MATMUL(MPtoC(:,:,ii,jj),WM(1:nVar,iGP,ii,jj))
              WP(1:nVar,iGP,ii,jj)=MATMUL(MPtoC(:,:,ii,jj),WP(1:nVar,iGP,ii,jj))
            END DO


          CASE DEFAULT
            ErrorMessage = "Error in ReconstructionFixY_INPUT_CONSERVED. ReconstructedVariable not available"
            WRITE(*,*) ErrorMessage
            STOP
        END SELECT





#ifdef PWCINTROUBLEDCELLS
#ifdef MULTIFLUID
        IF(PRESENT(Troubled_Cell_W)) THEN
          IF (Troubled_Cell_W(ii,jj) .NE. 0) THEN !*REPLACE BY PWC
            WM(1:nVar,nGPs,ii,jj) = W(1:nVar,ii,jj)
            WP(1:nVar,nGPs,ii,jj) = W(1:nVar,ii,jj)
          END IF
        END IF
#else
        PRINT*, "Impossible to use PWCINTROUBLEDCELLS"
        PRINT*, "Disable this flag"
        PRINT*, "You should not be here"
        STOP
#endif
#endif

#ifdef MULTIFLUID
        WM(nVar+1,1,ii,jj)=W(nVar,ii,jj)
        WP(nVar+1,1,ii,jj)=W(nVar,ii,jj)
#endif


      END DO
    END DO

  CASE(-1)

    PRINT*, "Reconstruction of Characteristic variables of Alex' and Shaoshuai not set up for subroutine ReconstructionFixY_INPUT_CONSERVED"
    STOP

  CASE DEFAULT
    ErrorMessage = "Error in ReconstructionFixY_INPUT. ReconstructedVariable not available"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT




!-------------------------------------------------------------------------------!
END SUBROUTINE ReconstructionFixY_INPUT_CONSERVED
!===============================================================================!
!
!
!
!===============================================================================!
FUNCTION MINMOD_3_INPUTS(x,y,z)
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: x, y, z
REAL            :: MINMOD_3_INPUTS
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL            :: absx, absy, absz
REAL            :: MA, MI
!-------------------------------------------------------------------------------!

  MINMOD_3_INPUTS = 0.0


  MA=MAX(x,y,z)
  MI=MIN(x,y,z)

  IF( MA*MI .LE. 0) THEN !*Different sign
    MINMOD_3_INPUTS = 0.0
  ELSE
    absx=ABS(x)
    absy=ABS(y)
    absz=ABS(z)

    !*SIGN(A,B) returns the value of A with the sign of B.
    !*NB: SIGN should work for reals
    MINMOD_3_INPUTS = SIGN(1.,x)*MIN(absx,absy,absz)

  END IF
!-------------------------------------------------------------------------------!
END FUNCTION MINMOD_3_INPUTS
!===============================================================================!
!
!
!
!===============================================================================!
FUNCTION MINMOD_4_INPUTS(w,x,y,z)
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: w, x, y, z
REAL            :: MINMOD_4_INPUTS
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL            :: absw, absx, absy, absz
REAL            :: MA, MI
!-------------------------------------------------------------------------------!

  MINMOD_4_INPUTS = 0.0


  MA=MAX(w,x,y,z)
  MI=MIN(w,x,y,z)

  IF( MA*MI .LE. 0) THEN !*Different sign
    MINMOD_4_INPUTS = 0.0
  ELSE
    absw=ABS(w)
    absx=ABS(x)
    absy=ABS(y)
    absz=ABS(z)

    !*SIGN(A,B) returns the value of A with the sign of B.
    !*NB: SIGN should work for reals
    MINMOD_4_INPUTS = SIGN(1.,x)*MIN(absw,absx,absy,absz)

  END IF
!-------------------------------------------------------------------------------!
END FUNCTION MINMOD_4_INPUTS
!===============================================================================!
!
!
!
!===============================================================================!
FUNCTION MINMOD_5_INPUTS(v,w,x,y,z)
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: v, w, x, y, z
REAL            :: MINMOD_5_INPUTS
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL            :: absv, absw, absx, absy, absz
REAL            :: MA, MI
!-------------------------------------------------------------------------------!

  MINMOD_5_INPUTS = 0.0


  MA=MAX(v,w,x,y,z)
  MI=MIN(v,w,x,y,z)

  IF( MA*MI .LE. 0) THEN !*Different sign
    MINMOD_5_INPUTS = 0.0
  ELSE
    absv=ABS(v)
    absw=ABS(w)
    absx=ABS(x)
    absy=ABS(y)
    absz=ABS(z)

    !*SIGN(A,B) returns the value of A with the sign of B.
    !*NB: SIGN should work for reals
    MINMOD_5_INPUTS = SIGN(1.,x)*MIN(absv,absw,absx,absy,absz)

  END IF
!-------------------------------------------------------------------------------!
END FUNCTION MINMOD_5_INPUTS
!===============================================================================!
!
!
!
!===============================================================================!
#ifdef ACTIVEFLUX
SUBROUTINE FromCellAveragesToPointValues()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
!*Staggering in X direction
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: W_qp_X

!*Staggering in Y direction
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: W_qp_Y
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

W_qp_X=0.
W_qp_Y=0.

!*Values for FLuxX
SELECT CASE (ABS(Reconstruction))
  CASE(1,10,2,20,21,22,23,24,25,26,27,28,29) !*Up to second order I can directly set point-values equal to cell averages. NB: :nGPs=1
    DO jj=1,nElemsY
      DO ii=0,nElemsX
        W_qp_X(1:nVar,1,ii,jj)=W_X(1:nVar,ii+1,jj) !*NB: There's a shift/staggering in the indices at the RHS
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Reconstruction not implemented in FromCellAveragesToPointValues"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!*Values for FLuxY
SELECT CASE (ABS(Reconstruction))
  CASE(1,10,2,20,21,22,23,24,25,26,27,28,29) !*Up to second order I can directly set point-values equal to cell averages. NB: :nGPs=1
    DO jj=0,nElemsY
      DO ii=1,nElemsX
        W_qp_Y(1:nVar,1,ii,jj)=W_Y(1:nVar,ii,jj+1) !*NB: There's a shift/staggering in the indices at the RHS
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Reconstruction not implemented in FromCellAveragesToPointValues"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT


!-------------------------------------------------------------------------------!
END SUBROUTINE FromCellAveragesToPointValues
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE AF_MUSCL(Q,QM,QP,WM,WP,dx,Qaverage)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)          :: dx
REAL,INTENT(IN)          :: Q(1:nVar,-1:1)
REAL,INTENT(IN)          :: QM(1:nVar)
REAL,INTENT(IN)          :: QP(1:nVar)
REAL,INTENT(OUT)         :: WM(1:nVar,1:nGPs)
REAL,INTENT(OUT)         :: WP(1:nVar,1:nGPs)
REAL,INTENT(IN),OPTIONAL :: Qaverage(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: sm, sp, sc, sc1, sc2, slope, sm1, sp1
INTEGER            :: iVar, iGP
REAL,    PARAMETER :: theta=1.0
REAL               :: W(1:nVar,-1:1), dx_input, slope_output(1:nVar)
!-------------------------------------------------------------------------------!

#if(1==0)

  W(1:nVar,0)=Q(1:nVar,0)

  !*Small stencil
  ! W(1:nVar,-1)=QM(1:nVar)
  ! W(1:nVar,+1)=QP(1:nVar)
  ! dx_input=dx/2.0

  !*Big stencil
  W(1:nVar,-1)=Q(1:nVar,-1)
  W(1:nVar,+1)=Q(1:nVar,+1)
  dx_input=dx


  CALL SBMMUSCL(&
          &  W(1:nVar,-1:1),&
          &  WM(1:nVar,1:nGPs),&
          &  WP(1:nVar,1:nGPs),&
          &  dx_input,slope_output=slope_output,&
          &  tau_input=0.5, theta_input=1.3)


#else


DO iGP=1,nGPs !*NB: nGPs is 1 in this case
  DO iVar = 1, nVar

#ifdef FOURINPUTSPOSTPROCESSING
    sm    = (QP(iVar) -  Q(iVar,0)) / (dx/2.0) * theta
    sc1   = (QP(iVar) - QM(iVar)) / dx
    sc2    = (Q(iVar,1) - Q(iVar,-1)) / (dx*2.0)
    sp    = ( Q(iVar,0) - QM(iVar)) / (dx/2.0) * theta

    slope = MINMOD_4_INPUTS(sm,sc1,sc2,sp)

#elif defined(ALINAPOSTPROCESSING)
    sp    = ( Qaverage(iVar) - QM(iVar)) / (dx/2.0) * theta    
    sc    = (QP(iVar) - QM(iVar)) / dx
    sm    = (QP(iVar) -  Qaverage(iVar)) / (dx/2.0) * theta
    slope = MINMOD_3_INPUTS(sm,sc,sp)
#else

    sp    = ( Q(iVar,0) - QM(iVar)) / (dx/2.0) * theta
#ifdef OTHERPOSTPROCESSING
    sc    = (Q(iVar,1) - Q(iVar,-1)) / (dx*2.0)
#else
    sc    = (QP(iVar) - QM(iVar)) / dx
#endif
    sm    = (QP(iVar) -  Q(iVar,0)) / (dx/2.0) * theta
    slope = MINMOD_3_INPUTS(sm,sc,sp)
#endif


    WM(iVar,iGP) = Q(iVar,0) - 0.5*slope*dx
    WP(iVar,iGP) = Q(iVar,0) + 0.5*slope*dx
  END DO
END DO

#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE AF_MUSCL
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE AF_MUSCL_Primitive_Based(Q,WM,WP,dx)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)          :: dx
REAL,INTENT(IN)          :: Q(1:nVar,-1:1)
REAL,INTENT(OUT)         :: WM(1:nVar,1:nGPs)
REAL,INTENT(OUT)         :: WP(1:nVar,1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: sm, sp, sc, sc1, sc2, slope, sm1, sp1
INTEGER            :: iVar, iGP
REAL,    PARAMETER :: theta=1.0
REAL               :: W(1:nVar,-1:1), dx_input, slope_output(1:nVar)
!-------------------------------------------------------------------------------!



DO iGP=1,nGPs !*NB: nGPs is 1 in this case
  DO iVar = 1, nVar

    sp    = ( Q(iVar,+1) - Q(iVar, 0) ) / dx * theta
    sc    = ( Q(iVar,+1) - Q(iVar,-1) ) / (dx*2.0)
    sm    = ( Q(iVar, 0) - Q(iVar,-1) ) / dx * theta

    slope = MINMOD_3_INPUTS(sm,sc,sp)


    WM(iVar,iGP) = Q(iVar,0) - 0.5*slope*dx
    WP(iVar,iGP) = Q(iVar,0) + 0.5*slope*dx
  END DO
END DO


!-------------------------------------------------------------------------------!
END SUBROUTINE AF_MUSCL_Primitive_Based
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE SIDED_AF_MUSCL(Q,QM,QP,WM,WP,dx,which_side)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)          :: dx
REAL,INTENT(IN)          :: Q(1:nVar,-1:1)
REAL,INTENT(IN)          :: QM(1:nVar)
REAL,INTENT(IN)          :: QP(1:nVar)
REAL,INTENT(OUT)         :: WM(1:nVar,1:nGPs)
REAL,INTENT(OUT)         :: WP(1:nVar,1:nGPs)
CHARACTER,INTENT(IN)     :: which_side
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: sm, sp, sc, sc1, sc2, slope, sm1, sp1
INTEGER            :: iVar, iGP
REAL,    PARAMETER :: theta=1.0
REAL               :: W(1:nVar,-1:1), dx_input, slope_output(1:nVar)
!-------------------------------------------------------------------------------!


DO iGP=1,nGPs !*NB: nGPs is 1 in this case
  DO iVar = 1, nVar


    sp    = ( Q(iVar,0) - QM(iVar)) / (dx/2.0) * theta
    sm    = (QP(iVar) -  Q(iVar,0)) / (dx/2.0) * theta

    IF (which_side .EQ. "C") THEN
      sc    = (Q(iVar,1) - Q(iVar,-1)) / (dx*2.0)
    ELSEIF ( which_side .EQ. "M" ) THEN
      sc    = (Q(iVar,0) - Q(iVar,-1)) / dx      
    ELSEIF ( which_side .EQ. "P" ) THEN
      sc    = (Q(iVar,+1) - Q(iVar,0)) / dx      
    ELSE
      sc=0.0
    END IF

    slope = MINMOD_3_INPUTS(sm,sc,sp)


    WM(iVar,iGP) = Q(iVar,0) - 0.5*slope*dx
    WP(iVar,iGP) = Q(iVar,0) + 0.5*slope*dx
  END DO
END DO


!-------------------------------------------------------------------------------!
END SUBROUTINE SIDED_AF_MUSCL
!===============================================================================!
!
!
!
!===============================================================================!
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
!===============================================================================!
!
!
!
!===============================================================================!
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
END MODULE MOD_Reconstruction
!-------------------------------------------------------------------------------!
