!===============================================================================!
MODULE MOD_JacobiIteration
!===============================================================================!

IMPLICIT NONE

PRIVATE

INTERFACE jacobi
  MODULE PROCEDURE jacobi
END INTERFACE

INTERFACE sparse_matmul
  MODULE PROCEDURE sparse_matmul
END INTERFACE

INTERFACE jacobi_2
  MODULE PROCEDURE jacobi_2
END INTERFACE


PUBLIC :: jacobi
PUBLIC :: sparse_matmul
PUBLIC :: jacobi_2

CONTAINS

 !*------------------------------------------------------------------------------------------------
 !* Compute the multiplication matrix * vector where matrix is given as CRS
 !*------------------------------------------------------------------------------------------------
 SUBROUTINE sparse_matmul(rowptr,colptr,matrixvalues, rhs, product)
   IMPLICIT NONE
   REAL, DIMENSION(:), INTENT(in)    :: matrixvalues, rhs
   INTEGER, DIMENSION(:), INTENT(in) :: rowptr, colptr
   REAL, DIMENSION(:), INTENT(inout) :: product

   INTEGER :: n, nnz, row, spIdx
 
   n=SIZE(rhs,DIM=1)
   nnz=SIZE(matrixvalues,DIM=1)
   product = 0.
   DO row=1,n
      DO spIdx=rowptr(row),rowptr(row+1)-1
         product(row)=product(row)+matrixvalues(spIdx)*rhs(colptr(spIdx))
      ENDDO
   ENDDO


  END SUBROUTINE sparse_matmul


  
 !*------------------------------------------------------------------------------------------------
 !* Jacobi iteration algorithm to find solution of (D+L)x=b
 !* The iterations are defined as x^(k+1)=D^-1(b-L*x^k)
 !* Input: L as a CRS (row, col, values), D as vector, b vector, sol vector output
 !*------------------------------------------------------------------------------------------------
 SUBROUTINE jacobi(rowptr,colptr,matrixvalues, diag, rhs, sol, iterations)
   IMPLICIT NONE
   REAL, DIMENSION(:), INTENT(in)    :: matrixvalues, rhs, diag
   INTEGER, DIMENSION(:), INTENT(in) :: rowptr, colptr
   REAL, DIMENSION(:), INTENT(inout) :: sol
   REAL, DIMENSION(SIZE(rhs,DIM=1))  :: solp
   REAL                              :: residual, tolerance
   INTEGER, INTENT(out)              :: iterations

   INTEGER :: n, nnz, k, maxIter, i
   n=SIZE(rhs,DIM=1)

   tolerance = 1.e-16
   maxIter = 1000
   
   k=0
   solp=rhs
   sol=0.
   residual = SUM(ABS(sol-solp))/n

   DO WHILE ( (tolerance<residual).AND. (k<maxIter))
      k=k+1
      CALL sparse_matmul(rowptr,colptr,matrixvalues,solp,sol)  !sol = L x^k
      DO  i=1, n
         sol(i)=(rhs(i)-sol(i))/diag(i)
      ENDDO 
      residual = SUM(ABS(sol-solp))/n
      !PRINT*, "res = ", residual

      solp=sol
   ENDDO

   IF (tolerance<residual) THEN
      print*, "Jacobi not converged. residual = ", residual
   ENDIF

   iterations = k

  END SUBROUTINE jacobi



 !*------------------------------------------------------------------------------------------------
 !* Jacobi iterative algorithm to find solution of Mx=b
 !* The iterations are defined as x^(k)=D^-1(b-N*x^(k-1))
 !* Input: M as a CRS (row, col, values), D as vector, b vector, sol vector output
 !*------------------------------------------------------------------------------------------------
 !*NB: Difference with respect to the pervious version is that here we give thw whole M in CRS
 !*------------------------------------------------------------------------------------------------
 SUBROUTINE jacobi_2(RowStart,Columns,Values,Diagonal,rhs,sol,iterations,initial_guess)
   IMPLICIT NONE
   !-------------------------------------------------------------------------------!
   ! >> FORMAL ARGUMENTS                                                           !
   !-------------------------------------------------------------------------------!
   INTEGER, DIMENSION(:), INTENT(IN)                         :: RowStart
   INTEGER, DIMENSION(:), INTENT(IN)                         :: Columns
   REAL,    DIMENSION(:), INTENT(IN)                         :: Values
   REAL,    DIMENSION(:), INTENT(IN)                         :: Diagonal
   REAL,    DIMENSION(:), INTENT(IN)                         :: rhs 
   REAL,    DIMENSION(:), INTENT(INOUT)                      :: sol
   INTEGER, INTENT(OUT)                                      :: iterations
   REAL,    DIMENSION(SIZE(sol,DIM=1)), INTENT(IN), OPTIONAL :: initial_guess
   !-------------------------------------------------------------------------------!
   ! >> LOCAL VARIABLES                                                            !
   !-------------------------------------------------------------------------------!
   INTEGER, PARAMETER                   :: maxIter = 1000000
   REAL,    PARAMETER                   :: tolerance = 1.e-12
   REAL,    DIMENSION(SIZE(sol,DIM=1))  :: solp
   REAL                                 :: residual
   INTEGER                              :: n
   INTEGER                              :: indr, indi, indc, indk
   REAL                                 :: drr
   REAL                                 :: sparseproduct
   LOGICAL, PARAMETER                   :: verbose = .FALSE.
   n=SIZE(sol,DIM=1)
   

   IF(PRESENT(initial_guess)) THEN
      solp=initial_guess
   ELSE
      solp=0.0
   END IF

   
   sol=0.0
   residual = SUM(ABS(sol-solp))/n

   indk=0

   DO WHILE ( (indk.EQ.0) .OR. ((tolerance<residual).AND. (indk<maxIter))) !*indk .EQ. 0 is essential to avoid glitches is initial_guess is made of zeros only

      indk=indk+1

      DO indr=1,n !*Loop on the rows
         drr=Diagonal(indr)
         sparseproduct=0.0
         DO indi=RowStart(indr),RowStart(indr+1)-1 !*Loop on the columns elements
            indc=Columns(indi) !*column index
            ! PRINT*, "row", indr, "index", indi, "column", indc
            IF (indc .NE. indr) THEN !*Not element on the diagonal
               sparseproduct=sparseproduct+Values(indi)*solp(indc) !*N*x^(k-1)
            END IF
         END DO
         sol(indr)=1.0/drr*(-sparseproduct+rhs(indr))
      END DO


      residual = SUM(ABS(sol-solp))/n
      IF ( verbose ) THEN
         PRINT*, "res = ", residual, "at iteration", indk
      END IF
      solp=sol

   ENDDO

   IF (tolerance<residual) THEN
      PRINT*, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      PRINT*, 'Jacobi did not converge'
      PRINT*, 'Iterations', maxIter
      PRINT*, 'Tolerance ', tolerance
      PRINT*, '  Final residual = ', residual
      PRINT*, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
   ENDIF

   iterations = indk


  END SUBROUTINE jacobi_2



END MODULE MOD_JacobiIteration
!-------------------------------------------------------------------------------!
