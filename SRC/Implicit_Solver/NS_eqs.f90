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
!> Routines to discretise and solve governing equations.
!> 1. Chebyshev pseudospectral method for  spatial discretization.
!> 2. The third-order backward difference method for temporal discretization.
!> 3. A fourth-order for-stage Runge-Kutta pseudo time integrator (RK) for solving the implicit discretized equations.
!=================================================================================================================================
MODULE MOD_NSEqs
USE MOD_Global
USE MOD_BC
IMPLICIT NONE
PRIVATE
SAVE
!---------------------------------------------------------------------------------------------------------------------------------
! PUBLIC VARIABLES
REAL(RP),PUBLIC:: Tao           !< Pseudo time step 
PUBLIC::Unsteady_Residual       !< Subroutine: calculate residual of the unsteady equations
PUBLIC::PseudoTimeIntegrator    !< Subroutine: pseudo time integrator 
!---------------------------------------------------------------------------------------------------------------------------------
! PRIVATE VARIABLES
INTEGER i,j,n
!=================================================================================================================================

CONTAINS
    
!=================================================================================================================================
!> Fourth-order for-stage Runge-Kutta pseudo time integrator (RK44) for solving the implicit discretized equations.
!=================================================================================================================================
SUBROUTINE PseudoTimeIntegrator
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(RP):: alpha(4)=(/0.25,1./3.,0.5,1.0/)   
REAL(RP):: uMax,vMax,LA,LB
INTEGER :: step
!=================================================================================================================================

! Compute pseudo time 
uMax=MAXVAL(ABS(Uone(:,:,2)))
vMax=MAXVAL(ABS(Uone(:,:,3)))
LA=(ABS(uMax)+SQRT(uMax**2+BetaP))/hx+1/(1./Kappa(2)*hx**2)
LB=(ABS(vMax)+SQRT(vMax**2+BetaP))/hy+1/(1./Kappa(2)*hy**2)
Tao=CFL/(LA+LB)
IF(Tao<1.D-14) Tao=0.1

! March for one pseudo time (RK44 solver)
Mid=Uone
DO step=1,4
    CALL Unsteady_Residual 
	DO n=1,nVar
		Mid(1:Nx-1,1:Ny-1,n)=Uone(1:Nx-1,1:Ny-1,n)+alpha(step)*Tao*Resi(1:Nx-1,1:Ny-1,n)
	ENDDO
	CALL BC_ADJ
ENDDO
Uone=Mid

END SUBROUTINE
    
!=================================================================================================================================
!> Calculate the residual of unsteady governing equaitons
!=================================================================================================================================
SUBROUTINE Unsteady_Residual
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(RP):: DU
!=================================================================================================================================

! Calculate the spatial terms
CALL SPACE_RESIDUAL

! Add temporal terms to the residual
DO n=2,nVar
		DO j=1,Ny-1; DO i=1,Nx-1  
			DU = 11.D0/6.D0*Mid(i,j,n)-3.D0*Uall(i,j,n,0) &
			      & +1.5D0*Uall(i,j,n,-1)-1.D0/3.D0*Uall(i,j,n,-2)
			Resi(i,j,n)=Resi(i,j,n)-DU/DT
		ENDDO; ENDDO  
ENDDO
Resi(0:Nx:Nx,0:Ny   ,:)=0
Resi(0:Nx   ,0:Ny:Ny,:)=0

END SUBROUTINE  
    
!=================================================================================================================================
!> Calculate the spatial terms
!=================================================================================================================================
SUBROUTINE SPACE_RESIDUAL
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(RP):: fxx(0:Nx),fyy(0:Ny)
REAL(RP):: alpha=1.D0,beta=0.D0
!=================================================================================================================================

! Calculate the convection terms 
CALL convect

! Add the diffusion terms to the residual
DO n=2,nVar
	    ! The second-order derivatives in x direction
		DO j=1,Ny-1
			!fxx
			CALL dgemv('n',Nx+1,Nx+1,alpha,DXX(0:Nx,0:Nx),Nx+1,Mid(0:Nx,j,n),1,beta,fxx(0:Nx),1)
			Resi(1:Nx-1,j,n)=Resi(1:Nx-1,j,n)+Kappa(n)*fxx(1:Nx-1)
		ENDDO
		
        ! The second-order derivatives in y direction
		DO i=1,Nx-1
			!fyy
			CALL dgemv('n',Ny+1,Ny+1,alpha,DYY(0:Ny,0:Ny),Ny+1,Mid(i,0:Ny,n),1,beta,fyy(0:Ny),1)
			Resi(i,1:Ny-1,n)=Resi(i,1:Ny-1,n)+Kappa(n)*fyy(1:Ny-1)
		ENDDO  
ENDDO

! Set the boundary residual to zero
Resi(0:Nx:Nx,0:Ny   ,:)=0
Resi(0:Nx   ,0:Ny:Ny,:)=0
END SUBROUTINE  

!=================================================================================================================================
!> Calculate the convection terms
!=================================================================================================================================
SUBROUTINE CONVECT
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(RP):: alpha=1.D0,beta=0.D0
REAL(RP):: phi_x(0:Nx,1:nVar) !dp_dx,du_dx,dv_dx,dT_dx
REAL(RP):: phi_y(0:Ny,1:nVar) !dp_dy,du_dy,dv_dy,dT_dy
!=================================================================================================================================

Resi=0.
! The first-order derivatives in x direction
DO j=1,Ny-1
	DO n=1,nVar
		CALL dgemv('n',Nx+1,Nx+1,alpha,DX(0:Nx,0:Nx),Nx+1,Mid(0:Nx,j,n),1,beta,phi_x(0:Nx,n),1)
	ENDDO
	Resi(1:Nx-1,j,1) = -BetaP*phi_x(1:Nx-1,2)
	Resi(1:Nx-1,j,2) = -Mid(1:Nx-1,j,2)*phi_x(1:Nx-1,2) - phi_x(1:Nx-1,1)
	Resi(1:Nx-1,j,3) = -Mid(1:Nx-1,j,2)*phi_x(1:Nx-1,3) + Mid(1:Nx-1,j,4)
	Resi(1:Nx-1,j,4) = -Mid(1:Nx-1,j,2)*phi_x(1:Nx-1,4)
ENDDO

! The first-order derivatives in x direction
DO i=1,Nx-1
    DO n=1,nVar
	    CALL dgemv('n',Ny+1,Ny+1,alpha,DY(0:Ny,0:Ny),Ny+1,Mid(i,0:Ny,n),1,beta,phi_y(0:Ny,n),1)
    ENDDO
	Resi(i,1:Ny-1,1)=Resi(i,1:Ny-1,1)-BetaP*phi_y(1:Ny-1,3)
	Resi(i,1:Ny-1,2)=Resi(i,1:Ny-1,2)-Mid(i,1:Ny-1,3)*phi_y(1:Ny-1,2)
	Resi(i,1:Ny-1,3)=Resi(i,1:Ny-1,3)-Mid(i,1:Ny-1,3)*phi_y(1:Ny-1,3)-phi_y(1:Ny-1,1)
	Resi(i,1:Ny-1,4)=Resi(i,1:Ny-1,4)-Mid(i,1:Ny-1,3)*phi_y(1:Ny-1,4)
ENDDO
! Set the boundary residual to zero
Resi(0:Nx:Nx,0:Ny   ,:)=0
Resi(0:Nx   ,0:Ny:Ny,:)=0
ENDSUBROUTINE

END MODULE
    
    
