!=================================================================================================================================
! Copyright (c) 2019  Mr. Wenqian Chen
! This file is part of PIDTS, a parallel-in-time implementation of the dual time stepping (DTS) method.
!
! PIDTS is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! PIDTS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. IF not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
!=================================================================================================================================
!> Routines to compute various errors, including:
!> 1. Convergence criterion--ErrorCriterion
!> 2. L1 error --ErrorL1
!> 3. L2 error --ErrorL2
!> 4. Max error--MaxError
!=================================================================================================================================
MODULE MOD_ErrorAnalysis
USE MOD_Global
IMPLICIT NONE
PUBLIC
!----------------------------------------------------------------------------------
! PUBLIC VARIABLES
REAL(RP):: ErrorCriterion
REAL(RP):: ErrorL1( nVar)
REAL(RP):: ErrorL2( nVar)
REAL(RP):: MaxError(nVar) 
!=================================================================================================================================

CONTAINS
!=================================================================================================================================
!> Compute convergence criterion
!=================================================================================================================================
SUBROUTINE Calerr
USE MOD_NSEqs, ONLY:Tao
IMPLICIT NONE
!----------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER i,j,n
!=================================================================================================================================
ErrorCriterion=0.
n=4
DO j=0,Ny; DO i=0,Nx
		ErrorCriterion=max(ErrorCriterion,abs(Uone(i,j,n)-U0(i,j,n)))
ENDDO;ENDDO
ErrorCriterion=ErrorCriterion/Tao
ErrorCriterion=log10( MAX(ErrorCriterion,1.E-30) )
ENDSUBROUTINE 

!=================================================================================================================================
!> Compute L1, L2 and Max errors 
!=================================================================================================================================
SUBROUTINE Calallerr
USE MOD_NSEqs
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER i,j,K,n
REAL(RP)::tem
!=================================================================================================================================
CALL Unsteady_Residual
ErrorL1(  1:nVar)=0.
ErrorL2(  1:nVar)=0.
MaxError( 1:nVar)=0.
DO n=1,nVar
	DO j=0,Ny; 	DO i=0,Nx
		ErrorL2(n) = ErrorL2(n) + Resi(i,j,n)**2.
		ErrorL1(n) = ErrorL1(n) + abs(Resi(i,j,n))
		IF(ABS(Resi(i,j,n)).GT.MaxError(n)) THEN
			MaxError(n)=ABS(Resi(i,j,n))
		END IF
	ENDDO; ENDDO
ENDDO
!PRINT*,'ErrorL2',ErrorL2(2:nVar)
ErrorL2=SQRT(ErrorL2/((Nx-1)*(Ny-1)))
ErrorL1=ErrorL1/((Nx-1)*(Ny-1))
ErrorL2=LOG10(ErrorL2)
ErrorL1=LOG10(ErrorL1)
END SUBROUTINE Calallerr
END MODULE