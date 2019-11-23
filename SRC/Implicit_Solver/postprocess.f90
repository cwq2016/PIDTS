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
!> Routines to compute the Nusselt number at left sidewall.
!=================================================================================================================================
MODULE MOD_Post
USE MOD_Global
USE MOD_Basis
IMPLICIT NONE
PRIVATE
SAVE
!---------------------------------------------------------------------------------------------------------------------------------
! PUBLIC VARIABLES
PUBLIC:: Init_Post      !< Routine to initialize Mod_post
PUBLIC:: ComputeNu      !< Routine to compute NuL
!---------------------------------------------------------------------------------------------------------------------------------
! PRIVATE VARIABLES
REAL(RP),ALLOCATABLE:: VDM_CGL2LG(:,:) !< Matrix to transform Chebyshev Gauss-Lobatto points
                                       !< to Legendre Gauss points
REAL(RP),ALLOCATABLE:: wLGP(:)         !< The integral coefficients of Legendre Gauss points
!=================================================================================================================================

CONTAINS

!=================================================================================================================================
!> Routines to compute VDM_CGL2LG and wLGP
!=================================================================================================================================
SUBROUTINE Init_Post
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(RP):: xCGP(0:Ny), wCGP(0:Ny), wCBary(0:Ny)
REAL(RP):: xLGP(0:Ny) 
!=================================================================================================================================
ALLOCATE( VDM_CGL2LG(0:Ny,0:Ny), wLGP(0:Ny))
CALL ChebyGaussLobNodesAndWeights(Ny,xCGP,wCGP)
xCGP(0:Ny)=xCGP(Ny:0:-1)
CALL BarycentricWeights(Ny,xCGP,wCBary)
CALL LegendreGaussNodesAndWeights(Ny,xLGP,wLGP)
CALL InitializeVandermonde(Ny,Ny,wCBary,xCGP,xLGP,VDM_CGL2LG)
END SUBROUTINE

!=================================================================================================================================
!> Routines to compute NuL
!=================================================================================================================================
FUNCTION ComputeNu() result(NuL)
USE MOD_Global
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER      :: i,j
REAL(RP)     :: NuL
REAL(RP)     :: NuR0(0:Ny),tmp(0:Ny),dthdx
!=================================================================================================================================

! Compute NuL at Chebyshev Gauss-Lobatto points
DO j=0,Ny
	NuR0(j)=-sum(DX(0,0:Nx)*Uone(0:Nx,j,4))
ENDDO

! Transform from Chebyshev Gauss-Lobatto points to Legendre Gauss points
tmp=0.d0
DO j=0,Ny;DO i=0,Ny
	tmp(i)=tmp(i)+VDM_CGL2LG(i,j)*NuR0(j)
ENDDO;ENDDO

! Integrate
NuL=sum(tmp*wLGP)*Leny/2.d0
END FUNCTION

END MODULE