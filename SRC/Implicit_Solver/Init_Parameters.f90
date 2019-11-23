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
!> The subroutine is used to read parameters from intput file, generate mesh and initialize parameters
!> for Chebyshev pseudospectral method 
!=================================================================================================================================
SUBROUTINE Init_Parameters
USE MOD_Global
USE MOD_BC
IMPLICIT NONE 
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER::i,j
!=================================================================================================================================

! Read parameters from input file
OPEN(100,FILE=TRIM(ADJUSTL(PATH_INPUT))//'/INPUT.txt',STATUS='OLD')

READ(100,NML=Rectangle)
READ(100,NML=Pulsating_Temperature)
READ(100,NML=Physical_Property)
READ(100,NML=Grid_Fine)
READ(100,NML=PTMU_Paras)
READ(100,NML=Unsteady_parameters)
READ(100,NML=Time_Parallel)
READ(100,NML=Output_Control)
    
Kappa=(/0.d0,sqrt(Pr/Ra),sqrt(Pr/Ra),1.d0/sqrt(Ra*Pr) /)

! Generate mesh 
xCoef=Lenx/2.0
yCoef=Leny/2.0
ALLOCATE(x(0:Nx),Y(0:Ny))
DO i=0,Nx
	x(i)=xCoef*cos(i*PI/Nx)
ENDDO
DO j=0,Ny
	y(j)=yCoef*cos(j*PI/Ny)
ENDDO


hx=x(0)-x(1)
hy=y(0)-y(1)

! Allocate space for solutions
ALLOCATE( Uone(0:Nx,0:Ny,nVar),U0(  0:Nx,0:Ny,nVar), &
		& Utmp(  0:Nx,0:Ny,nVar), &
		& Uall(0:Nx,0:Ny,nVar,-2:0), &
		& Mid(0:Nx,0:Ny,nVar),Resi(0:Nx,0:Ny,nVar))
U0  =0.D0
Uone=0.D0
Utmp=0.D0
Uall=0.D0
Mid=0.
Resi=0.

! Differntation matrix
ALLOCATE(DX(0:Nx,0:Nx),DXX(0:Nx,0:Nx))
ALLOCATE(DY(0:Ny,0:Ny),DYY(0:Ny,0:Ny))
CALL Velocity_Dev(Nx,DX,DXX,xCoef)   
CALL Velocity_Dev(Ny,DY,DYY,yCoef) 
	
CONTAINS


!=================================================================================================================================
!> Compute the first-order and second-order differentation matrix
!> (Sect. 3.3.4 -- Sect. 3.3.5, Peyret book)
!=================================================================================================================================
SUBROUTINE Velocity_Dev(n,vdev1,vdev2,coef)
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER ::n,i,j
REAL(RP)::vdev1(0:n,0:n),vdev2(0:n,0:n),coef
REAL(RP)::alpha,beta,sum
!=================================================================================================================================

alpha=1.0
beta=0.0

vdev1(0,0)=(2.*REAL(n)**2+1.0)/6.0
vdev1(n,n)=-(2.*REAL(n)**2+1.0)/6.0
DO i=0,n
	DO j=0,n
		IF(j.ne.i) then
			vdev1(i,j)=cl(i,n)*(-1)**(i+j)/&  
				(2*sin((j+i)*PI/(2*n))*sin((j-i)*PI/(2*n))*cl(j,n))
		endif
	ENDDO
	IF((i.ne.0).and.(i.ne.n)) then
		sum=0.
		DO j=0,n
			IF(j.ne.i) then
				sum=sum+vdev1(i,j)
			endif
		ENDDO
		vdev1(i,i)=-sum
	endif
ENDDO

CALL dgemm('n','n',n+1,n+1,n+1,alpha,vdev1,n+1,vdev1,n+1,beta,vdev2,n+1)

vdev1=vdev1/coef
vdev2=vdev2/(coef**2)

ENDSUBROUTINE
    
FUNCTION cl(i,n)
IMPLICIT NONE
INTEGER :: i,n
REAL(RP):: cl
IF((i.eq.0).or.(i.eq.n)) then
	cl=2.0
else
	cl=1.0
endif
RETURN
ENDFUNCTION
END SUBROUTINE Init_Parameters
    
    
