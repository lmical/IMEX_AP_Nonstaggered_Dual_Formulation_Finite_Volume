!===============================================================================!
MODULE MOD_Output
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE WriteMeshToDisk
  MODULE PROCEDURE WriteMeshToDisk
END INTERFACE

INTERFACE WriteSolutionToDisk
  MODULE PROCEDURE WriteSolutionToDisk
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: WriteMeshToDisk
PUBLIC :: WriteSolutionToDisk
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
SUBROUTINE WriteSolutionToDisk(OutputTime)
!-------------------------------------------------------------------------------!
USE MOD_equation,           ONLY: ExactFunction
USE MOD_equation,           ONLY: PrimToCons
USE MOD_equation,           ONLY: ConsToPrim
USE MOD_FiniteVolume2D_vars,ONLY: t
USE MOD_FiniteVolume2D_vars,ONLY: tGlobal
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Gravitational_Potential_Averages
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: ndims
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: VarNameVisu
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: WhichOutput
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGP
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
#if defined(CENTEREDPRIMITIVE)
USE MOD_FiniteVolume2D_vars,ONLY: Errors_PCSD
USE MOD_FiniteVolume2D     ,ONLY: Compute_Errors_PCSD
#endif
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE)
USE MOD_FiniteVolume2D_vars,ONLY: VarNameVisu_W
#endif
#ifdef CENTEREDPRIMITIVE
USE MOD_FiniteVolume2D_vars,ONLY: WC
#endif
#ifdef ACTIVEFLUX
USE MOD_FiniteVolume2D_vars,ONLY: W_X
USE MOD_FiniteVolume2D_vars,ONLY: W_Y
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary_X
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary_Y
USE MOD_Equation           ,ONLY: BoundaryConditions
USE MOD_Reconstruction     ,ONLY: First_Derivative_Central_Order2
#endif
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE) 
#ifdef IMEX
USE MOD_FiniteVolume2D_vars,ONLY: HyperbolicityLoss
#endif
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: OutputTime
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj, iGP, jGP, iVar
REAL,ALLOCATABLE   :: U_NVisu(:,:,:)
REAL,ALLOCATABLE   :: Error(:,:,:)
REAL,ALLOCATABLE   :: Uexact(:,:,:)
CHARACTER(LEN=255) :: FileNameTec,FileNameOct
REAL, DIMENSION(nVar, nGPs, nGPs) :: Utemp
REAL, DIMENSION(nVar) :: Conservation_Integral
REAL                  :: dxvx, dyvy
REAL                  :: int_abs_div_v_W_X, int_abs_div_v_W_Y
#ifdef CENTEREDPRIMITIVE
CHARACTER(LEN=255) :: FileNameOct_WC
CHARACTER(LEN=255) :: FileNameOct_WC_ERR
REAL,ALLOCATABLE   :: WC_NVisu(:,:,:)
REAL,ALLOCATABLE   :: WC_ERR_NVisu(:,:,:)
#endif
#ifdef ACTIVEFLUX
CHARACTER(LEN=255) :: FileNameOct_W_X
CHARACTER(LEN=255) :: FileNameOct_W_Y
REAL,ALLOCATABLE   :: W_X_NVisu(:,:,:)
REAL,ALLOCATABLE   :: W_Y_NVisu(:,:,:)
#endif
REAL, DIMENSION(nVar) :: WC_temp
REAL                  :: Error_temp(1:nVar), weight
INTEGER               :: indi, indj
!-------------------------------------------------------------------------------!

ALLOCATE(U_NVisu(1:nVar+1,1:nElemsX,1:nElemsY))
ALLOCATE(Error(1:nVar+1,1:nElemsX,1:nElemsY))
ALLOCATE(Uexact(1:nVar,1:nElemsX,1:nElemsY))
#ifdef CENTEREDPRIMITIVE
ALLOCATE(WC_NVisu(1:nVar+1,1:nElemsX,1:nElemsY))
ALLOCATE(WC_ERR_NVisu(1:nVar+1,1:nElemsX,1:nElemsY))
#endif
#ifdef ACTIVEFLUX
ALLOCATE(W_X_NVisu(1:nVar+1,1:nElemsX+1,1:nElemsY))
ALLOCATE(W_Y_NVisu(1:nVar+1,1:nElemsX,1:nElemsY+1))
#endif

FileNameOct = TIMESTAMP("SOLUTION_OCT",OutputTime)
FileNameOct = TRIM(FileNameOct)//".dat"
FileNameTec = TIMESTAMP("SOLUTION_TEC",OutputTime)
FileNameTec = TRIM(FileNameTec)//".dat"

#ifdef CENTEREDPRIMITIVE
FileNameOct_WC     = TRIM("WC_") // TRIM(FileNameOct)
FileNameOct_WC_ERR = TRIM("WC_Err_") // TRIM(FileNameOct)
#endif
#ifdef ACTIVEFLUX
FileNameOct_W_X = TRIM("W_X_") // TRIM(FileNameOct)
FileNameOct_W_Y = TRIM("W_Y_") // TRIM(FileNameOct)
#endif

#ifndef PRIMITIVEONLY
U_NVisu = 0.
#ifdef CENTEREDPRIMITIVE
WC_NVisu = 0.
WC_ERR_NVisu = 0.
#endif
Uexact = 0.
Utemp = 0.
Conservation_Integral=0.

#if defined(CENTEREDPRIMITIVE)
CALL Compute_Errors_PCSD()
#endif

DO jj=1,nElemsY
  DO ii=1,nElemsX
    U_NVisu(1:nVar,ii,jj) = U(1:nVar,ii,jj)
    U_NVisu(nVar+1,ii,jj) = Gravitational_Potential_Averages(ii,jj)
    DO iGP=1,nGPs
      DO jGP=1,nGPs
        CALL ExactFunction(InitialCondition, tGlobal, MeshGP(:,ii,jj,iGP,jGP), Utemp(1:nVar,iGP,jGP))
        Uexact(1:nVar, ii, jj) = Uexact(1:nVar, ii, jj)+ WeightsGP(iGP,jGP)* Utemp(1:nVar,iGP,jGP)
      END DO
    END DO
    Error(1:nVar,ii,jj) = U(1:nVar,ii,jj) - Uexact(1:nVar,ii,jj)
    Error(nVar+1,ii,jj) = 0. 
    Conservation_Integral=Conservation_Integral+U(1:nVar,ii,jj)
#ifdef CENTEREDPRIMITIVE
    WC_NVisu(1:nVar,ii,jj)     = WC(1:nVar,ii,jj)
    !*=================================================
    !*Error in Primitive variables
    ! CALL ConsToPrim(U(1:nVar,ii,jj),WC_temp(1:nVar))
    ! WC_ERR_NVisu(1:nVar,ii,jj) = ABS(WC(1:nVar,ii,jj)-WC_temp(1:nVar))
    !*=================================================
    !*Error in Conserved variables
    WC_ERR_NVisu(1:nVar,ii,jj) = Errors_PCSD(1:nVar,ii,jj)

#endif
  END DO
END DO
Conservation_Integral=Conservation_Integral*MESH_DX(1)*MESH_DX(2)
#endif

#ifdef ACTIVEFLUX
!*Computation of int_{\Omega} |div(v)| to check whether the div free constraints is satisfied
CALL BoundaryConditions(tGlobal)

!*Staggering in X direction
int_abs_div_v_W_X=0.0
DO jj=1,nElemsY
  DO ii=1,nElemsX+1
    dxvx=First_Derivative_Central_Order2(W_X(2,ii-1:ii+1,jj),MESH_DX(1))
    dyvy=First_Derivative_Central_Order2(W_X(3,ii,jj-1:jj+1),MESH_DX(2))
    int_abs_div_v_W_X=int_abs_div_v_W_X+ABS(dxvx+dyvy)
  END DO
END DO
int_abs_div_v_W_X=int_abs_div_v_W_X*MESH_DX(1)*MESH_DX(2)

!*Staggering in Y direction
int_abs_div_v_W_Y=0.0
DO jj=1,nElemsY+1
  DO ii=1,nElemsX
    dxvx=First_Derivative_Central_Order2(W_Y(2,ii-1:ii+1,jj),MESH_DX(1))
    dyvy=First_Derivative_Central_Order2(W_Y(3,ii,jj-1:jj+1),MESH_DX(2))
    int_abs_div_v_W_Y=int_abs_div_v_W_Y+ABS(dxvx+dyvy)
    ! PRINT*, ii, jj, dxvx+dyvy, W_Y(2,ii-1:ii+1,jj)
  END DO
END DO
int_abs_div_v_W_Y=int_abs_div_v_W_Y*MESH_DX(1)*MESH_DX(2)
#endif


PRINT*, "=========================="
PRINT*, "Time:", OutputTime 
#ifndef PRIMITIVEONLY
PRINT*, "Conservation integrals:", Conservation_Integral
PRINT*, "--------------------------"
DO iVar=1,nVar
  PRINT*, "Maximum and Minimum of conserved variable U", iVar, "is", MAXVAL(U(iVar,1:nElemsX,1:nElemsY)),  MINVAL(U(iVar,1:nElemsX,1:nElemsY))
END DO
#endif
#ifdef CENTEREDPRIMITIVE
PRINT*, "--------------------------"
DO iVar=1,nVar
  PRINT*, "Maximum and Minimum of conserved variable WC", iVar, "is", MAXVAL(WC(iVar,1:nElemsX,1:nElemsY)),  MINVAL(WC(iVar,1:nElemsX,1:nElemsY))
END DO
#endif
#ifdef ACTIVEFLUX
PRINT*, "--------------------------"
DO iVar=1,nVar
  PRINT*, "Maximum and Minimum of primitive variable W_X", iVar, "is", MAXVAL(W_X(iVar,1:nElemsX+1,1:nElemsY)),  MINVAL(W_X(iVar,1:nElemsX+1,1:nElemsY))
END DO
PRINT*, "--------------------------"
DO iVar=1,nVar
  PRINT*, "Maximum and Minimum of primitive variable W_Y", iVar, "is", MAXVAL(W_Y(iVar,1:nElemsX,1:nElemsY+1)),  MINVAL(W_Y(iVar,1:nElemsX+1,1:nElemsY+1))
END DO
PRINT*, "--------------------------"
PRINT*, "Divergence of velocity field in W_X", int_abs_div_v_W_X
PRINT*, "Divergence of velocity field in W_Y", int_abs_div_v_W_Y
PRINT*, "--------------------------"
#endif
#if defined(ACTIVEFLUX) || defined(CENTEREDPRIMITIVE) 
#ifdef IMEX
PRINT*, "Hyperbolicity loss", HyperbolicityLoss
#endif
#endif
PRINT*, "=========================="


#ifdef ACTIVEFLUX
DO jj=1,nElemsY
  DO ii=1,nElemsX+1 !*NB:+1 in X direction
    W_X_NVisu(1:nVar,ii,jj) = W_X(1:nVar,ii,jj)
    W_X_NVisu(nVar+1,ii,jj) = 0.0
#if MOMENTUMINPRIMITIVEVARIABLES
    W_X_NVisu(2:3,ii,jj) = W_X(2:3,ii,jj)/W_X(1,ii,jj)
#endif
  END DO
END DO


DO jj=1,nElemsY+1 !*NB:+1 in Y direction
  DO ii=1,nElemsX
    W_Y_NVisu(1:nVar,ii,jj) = W_Y(1:nVar,ii,jj)
    W_Y_NVisu(nVar+1,ii,jj) = 0.0
#if MOMENTUMINPRIMITIVEVARIABLES
    W_Y_NVisu(2:3,ii,jj) = W_Y(2:3,ii,jj)/W_Y(1,ii,jj)
#endif
  END DO
END DO
#endif


#ifdef CENTEREDPRIMITIVE
DO jj=1,nElemsY
  DO ii=1,nElemsX
    WC_NVisu(1:nVar,ii,jj) = WC(1:nVar,ii,jj)
    WC_NVisu(nVar+1,ii,jj) = 0.0
#if MOMENTUMINPRIMITIVEVARIABLES
    WC_NVisu(2:3,ii,jj) = WC(2:3,ii,jj)/WC(1,ii,jj)
#endif
  END DO
END DO
#endif



IF (WhichOutput .EQ. 1) THEN
#ifndef PRIMITIVEONLY
  CALL WriteDataToDisk_OCTAVE(&
        nVar,                &
        (/nElemsX,nElemsY/), &
        'XY',                &
        MeshBary(1:2,1:nElemsX,1:nElemsY),&
        VarNameVisu,         &
        U_NVisu,             &
        OutputTime,          &
        TRIM(FileNameOct),      &
        "FV2D")
#endif
#ifdef CENTEREDPRIMITIVE
  CALL WriteDataToDisk_OCTAVE(&
        nVar,                &
        (/nElemsX,nElemsY/), &
        'XY',                &
        MeshBary(1:2,1:nElemsX,1:nElemsY),&
        VarNameVisu_W,       &
        WC_NVisu,             &
        OutputTime,          &
        TRIM(FileNameOct_WC),      &
        "FV2D")

  CALL WriteDataToDisk_OCTAVE(&
        nVar,                &
        (/nElemsX,nElemsY/), &
        'XY',                &
        MeshBary(1:2,1:nElemsX,1:nElemsY),&
        VarNameVisu_W,         &
        WC_ERR_NVisu,             &
        OutputTime,          &
        TRIM(FileNameOct_WC_ERR),      &
        "FV2D")        
#endif
#ifdef ACTIVEFLUX
  CALL WriteDataToDisk_OCTAVE(                 &
        nVar,                                  &
        (/nElemsX+1,nElemsY/),                 & !*NB:+1 in X direction
        'XY',                                  &
        MeshBary_X(1:2,1:nElemsX+1,1:nElemsY), & !*NB:+1 in X direction
        VarNameVisu_W,                         &
        W_X_NVisu,        & !*NB:+1 in X direction
        OutputTime,                            &
        TRIM(FileNameOct_W_X),                 &
        "FV2D")

  CALL WriteDataToDisk_OCTAVE(                 &
        nVar,                                  &
        (/nElemsX,nElemsY+1/),                 & !*NB:+1 in Y direction
        'XY',                                  &
        MeshBary_Y(1:2,1:nElemsX,1:nElemsY+1), & !*NB:+1 in Y direction
        VarNameVisu_W,                         &
        W_Y_NVisu,        & !*NB:+1 in Y direction
        OutputTime,                            &
        TRIM(FileNameOct_W_Y),                 &
        "FV2D")        
#endif

ELSEIF (WhichOutput .EQ. 2) THEN
  CALL WriteDataToDisk_TECPLOT(&
        nVar,                &
        (/nElemsX,nElemsY/), &
        'XY',                &
        MeshBary(1:2,1:nElemsX,1:nElemsY),&
        VarNameVisu,         &
        U_NVisu,             &
        OutputTime,          &
        TRIM(FileNameTec),      &
        "FV2D")
ELSEIF (WhichOutput .EQ. 3) THEN
  CALL WriteDataToDisk_TECPLOT(&
        nVar,                &
        (/nElemsX,nElemsY/), &
        'XY',                &
        MeshBary(1:2,1:nElemsX,1:nElemsY),&
        VarNameVisu,         &
        U_NVisu,             &
        OutputTime,          &
        TRIM(FileNameTec),      &
        "FV2D")
    CALL WriteDataToDisk_OCTAVE(&
        nVar,                &
        (/nElemsX,nElemsY/), &
        'XY',                &
        MeshBary(1:2,1:nElemsX,1:nElemsY),&
        VarNameVisu,         &
        U_NVisu,             &
        OutputTime,          &
        TRIM(FileNameOct),      &
        "FV2D")
ENDIF


FileNameOct = TIMESTAMP("ERROR_OCT",OutputTime)
FileNameOct = TRIM(FileNameOct)//".dat"
FileNameTec = TIMESTAMP("ERROR_TEC",OutputTime)
FileNameTec = TRIM(FileNameTec)//".dat"
IF (WhichOutput .EQ. 1) THEN
  CALL WriteDataToDisk_OCTAVE(&
        nVar,                &
        (/nElemsX,nElemsY/), &
        'XY',                &
        MeshBary(1:2,1:nElemsX,1:nElemsY),&
        VarNameVisu,         &
        Error,               &
        OutputTime,          &
        TRIM(FileNameOct),      &
        "FV2D")
ELSEIF (WhichOutput .EQ. 2) THEN

  CALL WriteDataToDisk_TECPLOT(&
        nVar,                &
        (/nElemsX,nElemsY/), &
        'XY',                &
        MeshBary(1:2,1:nElemsX,1:nElemsY),&
        VarNameVisu,         &
        Error,               &
        OutputTime,          &
        TRIM(FileNameTec),      &
        "FV2D")
ELSEIF (WhichOutput .EQ. 3) THEN

  CALL WriteDataToDisk_TECPLOT(&
        nVar,                &
        (/nElemsX,nElemsY/), &
        'XY',                &
        MeshBary(1:2,1:nElemsX,1:nElemsY),&
        VarNameVisu,         &
        Error,               &
        OutputTime,          &
        TRIM(FileNameTec),      &
        "FV2D")
  CALL WriteDataToDisk_OCTAVE(&
        nVar,                &
        (/nElemsX,nElemsY/), &
        'XY',                &
        MeshBary(1:2,1:nElemsX,1:nElemsY),&
        VarNameVisu,         &
        Error,               &
        OutputTime,          &
        TRIM(FileNameOct),      &
        "FV2D")
END IF


DEALLOCATE(U_NVisu)
DEALLOCATE(Error)
DEALLOCATE(Uexact)
#ifdef CENTEREDPRIMITIVE
DEALLOCATE(WC_NVisu)
DEALLOCATE(WC_ERR_NVisu)
#endif
#ifdef ACTIVEFLUX
DEALLOCATE(W_X_NVisu)
DEALLOCATE(W_Y_NVisu)
#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE WriteSolutionToDisk
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WriteMeshToDisk()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: UNIT_FILE
#ifdef ACTIVEFLUX
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary_X
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary_Y
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
CHARACTER(LEN=255) :: FileName
!-------------------------------------------------------------------------------!

#ifndef PRIMITIVEONLY
Filename = TRIM("MESH_CoordinateX.dat")
WRITE(*,'(A22,1X,A)') &
  "Writing DATA to Disk: ", "OCTAVE ASCII"//" => "//TRIM(FileName)
OPEN(UNIT_FILE, FILE=Filename)
WRITE(UNIT_FILE,'(A)') "CoordinateX"
DO ii=1,nElemsX
  WRITE(UNIT_FILE,'(SP,1(ES21.14E2))') MeshBary(1,ii,1)
END DO
CLOSE(UNIT_FILE)

Filename = TRIM("MESH_CoordinateY.dat")
WRITE(*,'(A22,1X,A)') &
  "Writing DATA to Disk: ", "OCTAVE ASCII"//" => "//TRIM(FileName)
OPEN(UNIT_FILE, FILE=Filename)
WRITE(UNIT_FILE,'(A)') "CoordinateY"
DO jj=1,nElemsY
  WRITE(UNIT_FILE,'(SP,1(ES21.14E2))') MeshBary(2,1,jj)
END DO
CLOSE(UNIT_FILE)
#endif


#ifdef ACTIVEFLUX
Filename = TRIM("W_X_MESH_CoordinateX.dat")
WRITE(*,'(A22,1X,A)') &
  "Writing DATA to Disk: ", "OCTAVE ASCII"//" => "//TRIM(FileName)
OPEN(UNIT_FILE, FILE=Filename)
WRITE(UNIT_FILE,'(A)') "CoordinateX"
DO ii=1,nElemsX+1 !*NB:+1 in X direction
  WRITE(UNIT_FILE,'(SP,1(ES21.14E2))') MeshBary_X(1,ii,1)
END DO
CLOSE(UNIT_FILE)

Filename = TRIM("W_X_MESH_CoordinateY.dat")
WRITE(*,'(A22,1X,A)') &
  "Writing DATA to Disk: ", "OCTAVE ASCII"//" => "//TRIM(FileName)
OPEN(UNIT_FILE, FILE=Filename)
WRITE(UNIT_FILE,'(A)') "CoordinateY"
DO jj=1,nElemsY
  WRITE(UNIT_FILE,'(SP,1(ES21.14E2))') MeshBary_X(2,1,jj)
END DO
CLOSE(UNIT_FILE)

Filename = TRIM("W_Y_MESH_CoordinateX.dat")
WRITE(*,'(A22,1X,A)') &
  "Writing DATA to Disk: ", "OCTAVE ASCII"//" => "//TRIM(FileName)
OPEN(UNIT_FILE, FILE=Filename)
WRITE(UNIT_FILE,'(A)') "CoordinateX"
DO ii=1,nElemsX
  WRITE(UNIT_FILE,'(SP,1(ES21.14E2))') MeshBary_Y(1,ii,1)
END DO
CLOSE(UNIT_FILE)

Filename = TRIM("W_Y_MESH_CoordinateY.dat")
WRITE(*,'(A22,1X,A)') &
  "Writing DATA to Disk: ", "OCTAVE ASCII"//" => "//TRIM(FileName)
OPEN(UNIT_FILE, FILE=Filename)
WRITE(UNIT_FILE,'(A)') "CoordinateY"
DO jj=1,nElemsY+1 !*NB:+1 in Y direction
  WRITE(UNIT_FILE,'(SP,1(ES21.14E2))') MeshBary_Y(2,1,jj)
END DO
CLOSE(UNIT_FILE)
#endif


!-------------------------------------------------------------------------------!
END SUBROUTINE WriteMeshToDisk
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WriteDataToDisk_TECPLOT(nVar,nElems,CoordNames,Coords,VarNames,&
  OutputData,OutputTime,FileName,ZoneTitel)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: UNIT_FILE
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
INTEGER,INTENT(IN)          :: nVar
INTEGER,INTENT(IN)          :: nElems(2)
REAL,INTENT(IN)             :: OutputTime
REAL,INTENT(IN)             :: Coords(2,1:nElems(1),1:nElems(2))
REAL,INTENT(IN)             :: OutputData(1:nVar+1,1:nElems(1),1:nElems(2))
CHARACTER(LEN=2),INTENT(IN) :: CoordNames
CHARACTER(LEN=*),INTENT(IN) :: VarNames(1:nVar+1)
CHARACTER(LEN=*),INTENT(IN) :: FileName
CHARACTER(LEN=*),INTENT(IN) :: ZoneTitel
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
INTEGER                     :: ii, jj, iVar, Offset
INTEGER                     :: Width, Pos
CHARACTER(LEN=3)            :: space1
CHARACTER(LEN=7)            :: space2
CHARACTER(LEN=35)           :: VarString
CHARACTER(LEN=255)          :: FormatString
CHARACTER(LEN=512)          :: FormatTitle
CHARACTER(LEN=255)          :: temp1, temp2, temp3
REAL, DIMENSION(nVar+1)     :: tempOut
!-------------------------------------------------------------------------------!

WRITE(*,'(A22,1X,A)') &
  "Writing DATA to Disk: ", "TECPLOT ASCII"//" => "//TRIM(FileName)

Offset = 38
FormatTitle = ""
FormatTitle(1:37) = &
  'VARIABLES="Coordinate'//CoordNames(1:1)//'","Coordinate'//CoordNames(2:2)//'"'

DO iVar=1,nVar+1
  WRITE(VarString,'(A2,A,A1)') ',"', TRIM(VarNames(iVar)), '"'
  FormatTitle(Offset:Offset+LEN(TRIM(VarString))) = TRIM(VarString)
  Offset = Offset + LEN(TRIM(VarString))
END DO

WRITE(temp1,'(F15.8)') OutputTime
WRITE(temp2,'(I8)') nElems(1)
WRITE(temp3,'(I8)') nElems(2)
OPEN(UNIT_FILE, FILE = TRIM(FileName), STATUS = "REPLACE")
WRITE(UNIT_FILE,'(A)') FormatTitle(1:Offset-1)
WRITE(UNIT_FILE,'(A,A,A,A,A,A,A,A)') &
  'ZONE T="', TRIM(ZoneTitel), '",STRANDID=1,SOLUTIONTIME=', &
  TRIM(ADJUSTL(TRIM(temp1))), ', DATAPACKING=POINT, ZONETYPE=ORDERED, I=', &
  TRIM(ADJUSTL(TRIM(temp2))), ', J=', TRIM(ADJUSTL(TRIM(temp3)))

WRITE(FormatString,'(A4,I2,A15)') "(SP,", nVar+2+1, "(ES21.14E2,2X))"

DO jj=1,nElems(2)
  DO ii=1,nElems(1)
    DO iVar=1,nVar+1
      IF (ABS(OutputData(iVar,ii,jj))<1.e-70) THEN
        tempOut(iVar) = 0.d0
      ELSE
        tempOut(iVar) = OutputData(iVar,ii,jj)
      END IF
    ENDDO
    !WRITE(UNIT_FILE,FormatString) Coords(1:2,ii,jj), OutputData(1:nVar,ii,jj)
    WRITE(UNIT_FILE,FormatString) Coords(1:2,ii,jj), tempOut(1:nVar+1)
  END DO
END DO

CLOSE(UNIT_FILE)

!-------------------------------------------------------------------------------!
END SUBROUTINE WriteDataToDisk_TECPLOT
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WriteDataToDisk_OCTAVE(nVar,nElems,CoordNames,Coords,VarNames,&
  OutputData,OutputTime,FileName,ZoneTitel)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: UNIT_FILE
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
INTEGER,INTENT(IN)          :: nVar
INTEGER,INTENT(IN)          :: nElems(2)
REAL,INTENT(IN)             :: OutputTime
REAL,INTENT(IN)             :: Coords(2,1:nElems(1),1:nElems(2))
REAL,INTENT(IN)             :: OutputData(1:nVar+1,1:nElems(1),1:nElems(2))
CHARACTER(LEN=2),INTENT(IN) :: CoordNames
CHARACTER(LEN=*),INTENT(IN) :: VarNames(1:nVar+1)
CHARACTER(LEN=*),INTENT(IN) :: FileName
CHARACTER(LEN=*),INTENT(IN) :: ZoneTitel
!-------------------------------------------------------------------------------!
! LOCAL VARIABLES
!-------------------------------------------------------------------------------!
INTEGER                     :: ii, jj, iVar, Offset
INTEGER                     :: Width, Pos
CHARACTER(LEN=3)            :: space1
CHARACTER(LEN=7)            :: space2
CHARACTER(LEN=35)           :: VarString
CHARACTER(LEN=255)          :: FormatString
CHARACTER(LEN=512)          :: FormatTitle
CHARACTER(LEN=255)          :: temp1, temp2, temp3
REAL, DIMENSION(nVar+1)       :: tempOut
!-------------------------------------------------------------------------------!

WRITE(*,'(A22,1X,A)') &
  "Writing DATA to Disk: ", "OCTAVE ASCII"//" => "//TRIM(FileName)

Pos   = 1
Width = 23
space1 = "   "
FormatTitle = ""

DO iVar=1,nVar+1
  WRITE(VarString,'(A)') TRIM(VarNames(iVar))
  FormatTitle(Pos:Pos+Width-1) = ADJUSTL(TRIM(VarString))
  Pos = Pos + Width
END DO

OPEN(UNIT_FILE, FILE = TRIM(FileName), STATUS = "REPLACE")
WRITE(UNIT_FILE,'(A)') FormatTitle(1:Pos-1)
WRITE(FormatString,'(A4,I2,A15)') "(SP,", nVar+1, "(ES21.14E2,2X))"

DO jj=1,nElems(2)
  DO ii=1,nElems(1)
    DO iVar=1,nVar+1
      IF (ABS(OutputData(iVar,ii,jj))<1.e-70) THEN
        tempOut(iVar) = 0.d0
      ELSE
        tempOut(iVar) = OutputData(iVar,ii,jj)
      END IF
    ENDDO
    WRITE(UNIT_FILE,FormatString) tempOut
  END DO
END DO

CLOSE(UNIT_FILE)

!-------------------------------------------------------------------------------!
END SUBROUTINE WriteDataToDisk_OCTAVE
!===============================================================================!
!
!
!
!===============================================================================!
FUNCTION TIMESTAMP(Filename,Time)
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,            INTENT(IN) :: Time
CHARACTER(LEN=*),INTENT(IN) :: Filename
CHARACTER(LEN=255)          :: TimeStamp
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER                     :: i
!-------------------------------------------------------------------------------!

WRITE(TimeStamp,'(F15.7)') Time
DO i=1,LEN(TRIM(TimeStamp))
  IF(TimeStamp(i:i) .EQ. " ") TimeStamp(i:i) = "0"
END DO
TimeStamp=TRIM(Filename)//'_'//TRIM(TimeStamp)

!-------------------------------------------------------------------------------!
END FUNCTION TIMESTAMP
!===============================================================================!
!
!
!
!===============================================================================!
END MODULE MOD_Output
!-------------------------------------------------------------------------------!
